source("utils.R")

source("wolfel.R")
source("he.R")
source('kucirka_fitting.R')

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  mutate(syndromic_sensitivity = syndromic_sensitivity)  %>%
  bind_cols(., list(
    `low` = 
      crossing(screening = c(TRUE,FALSE),
               first_test_delay = 0,
               second_test_delay = NA), 
    `moderate` = 
      crossing(screening = c(TRUE,FALSE),
               first_test_delay = c(3,5,7),
               second_test_delay = NA),
    `high` = 
      crossing(screening = TRUE,
               first_test_delay =  c(0:3),
               second_test_delay = c(2,4,6)),
    `maximum` = 
      crossing(screening = c(TRUE,FALSE),
               first_test_delay = 14,
               second_test_delay = NA)) %>%
      bind_rows(.id = "stringency")) %>% 
  mutate(scenario=row_number()) 

notification_t <- seq(from=0,to=12,by=2)

run_analysis <- 
  function(n_arrival_sims  = 1000,
           n_sec_cases     = 10,
           seed            = 145,
           notification_t){
    
    browser()
    set.seed(seed)
    
    #Parameters
    
    incubation_times <- make_incubation_times(
      n_travellers = n_sec_cases,
      pathogen     = pathogen,
      syndromic_sensitivity = unique(input$syndromic_sensitivity))
    
    #exposure date relative to index cases exposure
    incubation_times %<>% mutate(exposed_t=si$r(n()))
    
    #cross with scenarios
    incubation_times %<>% crossing(input) 
    
    #cross with notification times then remove exposures post-notification
    incubation_times %<>% 
      crossing(notification_t) %>% 
      mutate(remove=ifelse(exposed_t>notified_t,TRUE,FALSE)) %>% 
      dplyr::filter(!remove) 
  
    #calc_outcomes 
    incubation_times %<>% calc_outcomes(.,dat_gam)
    
    #when released
    incubation_times %<>% when_released()
    
    return(incubation_times)
    
  }


# 1000 secondary cases with draw from serial interval. 
# 0 t is time of exposure of first case
sec_cases <- tibble(i=1:1000) %>% 
            mutate(exposed_t=si$r(i))

#
sec_cases %<>% 
  #cross with notification times 
  crossing(notified_t=c(2,4,6,8,10,12)) %>% 
  #cross with scenarios
  crossing(input) %>% 
  #remove exposures occurring after notification
  mutate(remove=ifelse(exposed_t>notified_t,TRUE,FALSE)) %>% 
  dplyr::filter(!remove) %>% 
  #calc incubation times
  make_incubation_times()
  #calc outcomes
  calc_outcomes(.,dat_gam)



