# parameters for simulation

# Ashcroft et al. infectivity profile
# https://smw.ch/article/doi/smw.2020.20336
infect_shape = 97.18750 
infect_rate  =  3.71875
infect_shift = 25.62500

# Lauer et al. incubation period
# https://www.acpjournals.org/doi/10.7326/M20-0504
rriskDistributions::get.lnorm.par(q = c(5.1, 11.5),
                                  p = c(0.5, 0.975),
                                  plot = F,
                                  show.output = F) %>%
  as.list %>%
  set_names(., c("mu_inc", "sigma_inc")) -> inc_parms

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      # variance calculated by inverting confidence interval
      list(mu_inf    = 9.1,
           sigma_inf = 14.7)),
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

waning_none <- function(x){
  waning_points(x, X = 0, Y = 1)
}

adhere_10 <- function(x){
  waning_points(x, X = 0, Y = 0.1)
}
adhere_20 <- function(x){
  waning_points(x, X = 0, Y = 0.2)
}
adhere_30 <- function(x){
  waning_points(x, X = 0, Y = 0.3)
}
adhere_40 <- function(x){
  waning_points(x, X = 0, Y = 0.4)
}
adhere_50 <- function(x){
  waning_points(x, X = 0, Y = 0.5)
}
adhere_60 <- function(x){
  waning_points(x, X = 0, Y = 0.6)
}
adhere_70 <- function(x){
  waning_points(x, X = 0, Y = 0.7)
}
adhere_80 <- function(x){
  waning_points(x, X = 0, Y = 0.8)
}
adhere_90 <- function(x){
  waning_points(x, X = 0, Y = 0.9)
}
adhere_100 <- function(x){
  waning_points(x, X = 0, Y = 1)
}

waning_constant <- function(x){
  waning_points(x, X = 0, Y = 0.75)
}

# waning_drop <- function(x){
#   waning_piecewise_linear(x, 0.75, 0.25, 7, 14)
# }

# waning_linear <- function(x){
#   waning_piecewise_linear(x, ymax = 0.75, .16, 0, 8.3)
# }

# waning_canada_community <- function(x){
#   waning_points(x, X = c(0, 30), Y = c(1, 0.541), log = T)
# }

waning_canada_total <- function(x){
  waning_points(x, X = c(0, 30), Y = c(1, 0.158), log = T)
}

smith_uk <- function(x){
  waning_points(x, X = 0, Y = 0.109)
}

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `No testing` = 
      data.frame(sampling_freq=NA,
               tests=F,
               multiple_tests=NA,
               assay=NA,
               quar_dur=10),
    `Daily testing` = 
      data.frame(sampling_freq=1,
               tests=T,
               multiple_tests=T,
               assay="LFA",
               quar_dur=NA), 
    `Test at end of quarantine only` = 
      data.frame(sampling_freq=NA,
               multiple_tests=F,
               assay="LFA",
               quar_dur = 10)) %>% 
      bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window =  10,
           index_test_delay    =  c(1, 2, 3),  # time to entering quarantine
           delay_scaling       =  c(1, 0.5, 0),
           waning              =c("adhere_100"),
           adherence=c(0,0.5,1)) %>% 
  mutate(test_to_tracing=4*delay_scaling) %>% 
  mutate(scenario=row_number()) %>% 
  #calculate time until release from exposure for each scenario
  mutate(time_since_exp=quar_dur)

curves <- readr::read_csv(here::here("data","curve_out.csv"))

curves %<>% select(idx = iter,
                   days_since_infection = diff,
                   value)

curves %>% mutate(days_since_infection=days_since_infection*2/3) %>% sample_n(1000) %>% ggplot()+geom_point(aes(x=days_since_infection,y=value,group=idx),alpha=0.2)
