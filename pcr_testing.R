library(here)
setwd(here())

set.seed(145)
source('packages.R')

if (Sys.info()["nodename"] == "Sams-MacBook-Pro.local"){
  future::plan(multicore)
} else {
  future::plan(multiprocess)
}

source('utils.R')

#Parameters
n_arrival_sims <- 1000

tic()

# 
# input %<>% filter(stringency == "low" |
#                     stringency == "maximum" |
#                     (stringency %in% c("moderate", "high") &
#                        first_test_delay == 3 & (pre_board_screening != 4 |
#                                                   is.na(pre_board_screening)) ) )

n_travellers <- 1e4

incubation_times <- crossing(idx  = 1:n_travellers,
                             type = c("symptomatic",
                                      "asymptomatic") %>%
                               factor(x = .,
                                      levels = .,
                                      ordered = T)) %>%
  split(.$type) %>%
  map2_df(.x = .,
          .y = pathogen,
          ~mutate(.x,
                  exp_to_onset   = time_to_event(n = n(),
                                                 mean = .y$mu_inc, 
                                                 var  = .y$sigma_inc),
                  onset_to_recov = time_to_event(n = n(),
                                                 mean = .y$mu_inf, 
                                                 var  = .y$sigma_inf))) 

source("wolfel.R")
source("he.R")
# infectious period from time of onset to no longer infectious
incubation_times %<>% 
  mutate(u = runif(n = nrow(.), 0.01, 0.99)) %>%
  mutate(inf_from_onset = 
           approx(x    = wolfel_pred$y, 
                  y    = wolfel_pred$day, 
                  xout = u)$y,
         pre_symp_lead  = 
           approx(x    = HE$p,
                  y    = HE$delay,
                  xout = pmin(1 - 1e-5,
                              pmax(1e-5,
                                   pgamma(q = exp_to_onset,
                                          shape = inc_parms$shape,
                                          scale = inc_parms$scale))))$y
  )

incubation_times %<>% 
  mutate(onset     = exp_to_onset,
         inf_start = onset - pre_symp_lead,
         inf_end   = ifelse(type == "asymptomatic",
                            exp_to_onset + onset_to_recov,
                            exp_to_onset + inf_from_onset),
         symp_end  = ifelse(type == "asymptomatic",
                            onset, # but really never matters because asymptomatics are never symptomatic!
                            exp_to_onset + onset_to_recov),
         inf_dur   = inf_end - inf_start,
         symp_dur  = symp_end - onset)

# add flight
incubation_times %<>% 
  mutate(flight_departure = runif(n = nrow(.),
                                  min = 0,
                                  max = onset + onset_to_recov)) %>%
  mutate(symp_screen_label      = 
           flight_departure > onset &
           flight_departure < symp_end &
           runif(n = nrow(.)) < unique(input$syndromic_sensitivity))

incubation_times %<>% gen_screening_draws

# Infected arrivals

inf_arrivals <- as.list(c("EU", "USA")) %>%
  set_names(., .) %>%
  map(~make_prevalence(prev_est_region = prev_est_region,
                       origin_country = .x, 
                       n = n_arrival_sims)) %>%
  map_dfr(.id = "country", ~data.frame(pi = .x)) %>%
  mutate(alpha = rbeta(n = nrow(.),
                       shape1 = asymp_fraction$shape1,
                       shape2 = asymp_fraction$shape2)) %>%
  nest(data = -c(country)) %>%
  inner_join(dplyr::filter(flight_vols, year == 2020) %>%
               gather(country, trav_vol),
             by = "country") %>%
  mutate(trav_vol = 10000)


inf_arrivals <-  
  mutate(inf_arrivals,
         travellers = map2(.x = data,
                           .y = trav_vol,
                           .f = 
                             ~make_travellers(
                               x = .x,
                               incubation_times = incubation_times,
                               trav_vol = .y,
                               trav_vol_p = 1,
                               fixed = TRUE))) 


inf_arrivals  <- 
  mutate(inf_arrivals,
         individuals = future_map(.x = travellers,
                                  .f = travellers_to_individuals,
                                  incubation_times = incubation_times))

inf_arrivals %<>% unnest(individuals)

# add flight duration
inf_arrivals %<>% inner_join(flight_times)


## table 1 flight vol data

inf_arrivals %>%
  mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
  count(country, type, sim, .drop = F) %>%
  group_by_at(.vars = vars(-n, -sim)) %>%
  nest() %>%
  mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
  unnest_wider(Q) %>%
  select(-data) %>%
  group_by(country, type) %>%
  transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
  mutate(type = str_to_title(type)) %>%
  spread(country, value)

## totals for main text
inf_arrivals %>%
  mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
  count(country, sim, .drop = F) %>%
  group_by_at(.vars = vars(-n, -sim)) %>%
  nest() %>%
  mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
  unnest_wider(Q) %>%
  select(-data) %>%
  group_by(country) %>%
  transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
  spread(country, value)

source('kucirka_fitting.R')

# Cross arrivals with scenarios
arrival_scenarios <- crossing(input, inf_arrivals)

# calculate outcomes of screening
arrival_scenarios %<>% calc_outcomes(., dat_gam)

# Calculate when released
arrival_released <- when_released(arrival_scenarios)


# Calculate stage of infectiousness when released
arrival_released_times <- arrival_released %>%
  # select(-symp_screen_label) %>%
  # tidyr::unnest(released) %>% 
  # select(-c(first_test_delay,second_test_delay)) %>% 
  # unnest(input) %>% 
  mutate(stage_released = as.factor(case_when(
    is.na(released_t)         ~ "Prevented from boarding",
    released_t < inf_end      ~ "Infectious",
    released_t >= inf_end     ~ "Post-infectious"
  ))) %>% 
  mutate(days_released_inf = inf_end-released_t)

toc()
# 
# saveRDS(arrival_released_times,
# 
#         here("data",
#              paste0(paste("arrival",
#                           "released",
#                           "times",
#                           "10000",
#                           "trav",
#                           n_arrival_sims,
#                           "sims",
#                           sub(x = gsub(x = Sys.time(),
#                                        pattern = "(-|:)",
#                                        replacement = ""),
#                               pattern = "\\s",
#                               replacement = "_"),
#                           sep = "_"),".RDS"))
# )

