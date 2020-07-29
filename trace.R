source("utils.R")

source("wolfel.R")
source("he.R")

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

notification_times <- tibble(notified_t=seq(from=0,to=12,by=2))

run_analysis <- 
  function(n_arrival_sims  = 1000,
           n_sec_cases     = 10,
           seed            = 145,
           notification_t = notification_times){
    
    #browser()
    set.seed(seed)
    
    #Parameters
    
    incubation_times <- make_incubation_times(
      n_travellers = n_sec_cases,
      pathogen     = pathogen,
      syndromic_sensitivity = unique(input$syndromic_sensitivity))
    
    #exposure date relative to index cases exposure
    #BEAR IN MIND NEED TO SHIFT ALL t BY EXPOSURE DATE
    incubation_times %<>% 
      mutate(exposed_t= si$r(n())) %>% 
      #obviously some better way to do this
      mutate(onset    = onset+exposed_t,
             inf_start= inf_start+exposed_t,
             inf_end  = inf_end+exposed_t,
             symp_end = symp_end+exposed_t)
    

    
    #cross with scenarios
    incubation_times %<>% crossing(input) 
    
    #cross with notification times then remove exposures post-notification
    incubation_times %<>% 
      crossing(notification_t) %>% 
      mutate(remove=ifelse(exposed_t>notified_t,TRUE,FALSE)) %>% 
      dplyr::filter(!remove) 
    
    source('kucirka_fitting.R',local=T)  
    
    #calc outcomes 
    incubation_times %<>% calc_outcomes(.,dat_gam)
    
    #when released
    incubation_times %<>% when_released()
    
    return(incubation_times)
    
  }

results <- run_analysis()

results %>% 
  


