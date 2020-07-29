source("packages.R")

if (Sys.info()["nodename"] == "Sams-MacBook-Pro.local"){
  future::plan(multicore)
} else {
  future::plan(multiprocess)
  options(future.globals.maxSize = 2000*1024^2 )
}

source("utils.R")

# should really ask for an incubation_times, which is common across all
# should really ask for input, common across all
run_analysis <- 
  function(n_arrival_sims  = 1000,
           countries       = c("EU", "USA"),
           trav_vol_manual = NULL,
           n_travellers    = 1e4,
           trav_vol_p      = 1,
           fixed           = TRUE,
           seed            = 145){
    
    #browser()
    set.seed(seed)
    #Parameters
    
    incubation_times <- make_incubation_times(
      n_travellers = n_travellers,
      pathogen     = pathogen,
      syndromic_sensitivity = unique(input$syndromic_sensitivity))
    
    # Infected arrivals
    
    inf_arrivals <- make_inf_arrivals(
      countries       = countries,
      prev_est_region = prev_est_region,
      n_arrival_sims  = n_arrival_sims,
      asymp_fraction  = asymp_fraction,
      trav_vol_p      = trav_vol_p,
      flight_vols     = flight_vols,
      flight_times    = flight_times,
      trav_vol_manual = trav_vol_manual,
      incubation_times = incubation_times,
      fixed            = fixed)
    
    
    # Cross arrivals with scenarios
    arrival_scenarios <- make_arrival_scenarios(input, inf_arrivals, incubation_times)
    
    # Calculate when released
    arrival_released <- when_released(arrival_scenarios)
    
    # Calculate stage of infectiousness when released
    arrival_released_times <- stage_when_released(arrival_released)
    
    
    
    return(arrival_released_times)
    
  }

run_real_trav_vol_analysis <- function(n_arrival_sims = 1000,
                                       trav_vol_p=7/30,
                                       log_scale = FALSE,
                                       fixed = FALSE,
                                       text_size=2.5,
                                       ylabA = "Number of infectious travellers released per week",
                                       ylabB = "Remaining infectivity (person-days)\nof week's released travellers"){
  
  #browser()
  set.seed(145)
  
  #Parameters
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
    mutate(trav_vol = round(trav_vol/2))
  
  
  inf_arrivals <-  
    mutate(inf_arrivals,
           travellers = map2(.x = data,
                             .y = trav_vol,
                             .f = 
                               ~make_travellers(
                                 x = .x,
                                 incubation_times = incubation_times,
                                 trav_vol = .y,
                                 trav_vol_p = trav_vol_p,
                                 fixed = fixed))) 
  
  inf_arrivals  <- 
    mutate(inf_arrivals,
           individuals = future_map(.x = travellers,
                                    .f = travellers_to_individuals,
                                    incubation_times = incubation_times))
  
  inf_arrivals %<>% unnest(individuals)
  
  # add flight duration
  inf_arrivals %<>% inner_join(flight_times)
  
  ## table 1 flight vol data
  
  table_1 <- inf_arrivals %>%
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
  main_text_dat <- inf_arrivals %>%
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
  
  source('kucirka_fitting.R',local = T)
  
  # Cross arrivals with scenarios
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  # Calculate when released
  arrival_released <- when_released(arrival_scenarios)
  
  # Calculate stage of infectiousness when released
  arrival_released_times <- stage_when_released(arrival_released)
  
  
  source("make_plots_alt.R",local = T)
  
  list("png", "pdf") %>%
    map(~ggsave(filename = paste0("results/figS2_v2.",.x),
                plot = fig3,
                width = 210, height = 280, units="mm",
                dpi = 320))
  
  return(list(table_1 = table_1,
              main_text_dat = main_text_dat,
              n_data = fig3A_data,
              pd_data = fig3B_data,
              arrival_released_times = arrival_released_times))
  
}


run_rr_analysis <- function(n_arrival_sims = 500,
                            trav_vol_manual=100000, 
                            baseline_scenario = "7",
                            text_size = 2.5,
                            log_scale = FALSE,
                            fixed = TRUE){
  set.seed(145)
  browser()
  #Parameters
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
    mutate(trav_vol = trav_vol_manual)
  
  
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
                                 fixed = fixed))) 
  
  inf_arrivals  <- 
    mutate(inf_arrivals,
           individuals = future_map(.x = travellers,
                                    .f = travellers_to_individuals,
                                    incubation_times = incubation_times))
  
  inf_arrivals %<>% unnest(individuals)
  
  # add flight duration
  inf_arrivals %<>% inner_join(flight_times)
  
  source('kucirka_fitting.R',local = T)
  
  gc()
  # Cross arrivals with scenarios
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  # Calculate when released
  arrival_released <- when_released(arrival_scenarios)
  
  # Calculate stage of infectiousness when released
  arrival_released_times <- stage_when_released(arrival_released)
  
  
  source("rr_calc.R",local = T)
  
  list("png", "pdf") %>%
    map(~ggsave(filename = paste0("results/baseline_",baseline_scenario,"_rr_figs.",.x),
                plot=rr_figs,
                width = 210, height = 280, units="mm",
                dpi = 320))
  
  write.csv(n_fig_data,paste0("results/baseline_",baseline_scenario,"_n_RR_results.csv"))
  write.csv(pd_fig_data,paste0("results/baseline_",baseline_scenario,"_pd_RR_results.csv"))
  
  return(list(arrival_released_times = arrival_released_times,n_fig_data=n_fig_data,pd_fig_data=pd_fig_data))
}

run_cumul_pd_analysis <- 
  function(n_sims=1000,
           trav_vol_manual=10000,
           trav_vol_p=1,
           baseline_scenario = "7",
           log_scale = FALSE,
           text_size = 2.5,
           fixed = TRUE, 
           ylabA = "Cumulative remaining infectivity \n(person-days) of released travellers"){
    
    browser()
    set.seed(145)
    
    #Parameters
    n_arrival_sims <- n_sims
    
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
      mutate(trav_vol = trav_vol_manual)
    
    
    inf_arrivals <-  
      mutate(inf_arrivals,
             travellers = map2(.x = data,
                               .y = trav_vol,
                               .f = 
                                 ~make_travellers(
                                   x = .x,
                                   incubation_times = incubation_times,
                                   trav_vol = .y,
                                   trav_vol_p = trav_vol_p,
                                   fixed = fixed))) 
    
    inf_arrivals  <- 
      mutate(inf_arrivals,
             individuals = future_map(.x = travellers,
                                      .f = travellers_to_individuals,
                                      incubation_times = incubation_times))
    
    inf_arrivals %<>% unnest(individuals)
    
    # add flight duration
    inf_arrivals %<>% inner_join(flight_times)
    
    source('kucirka_fitting.R',local=T)
    
    # Cross arrivals with scenarios
    arrival_scenarios <- crossing(input, inf_arrivals)
    
    # calculate outcomes of screening
    arrival_scenarios %<>% calc_outcomes(., dat_gam)
    
    # Calculate when released
    arrival_released <- when_released(arrival_scenarios)
    
    
    # Calculate stage of infectiousness when released
    arrival_released_times <- stage_when_released(arrival_released)
    
    
    source("make_plots_pd.R",local = T)
    
    list("png", "pdf") %>%
      map(~ggsave(filename = paste0("results/figS5_v2.",.x),
                  plot = figS5,
                  width = 210, height = 280, units="mm",
                  dpi = 200))
    
  }
