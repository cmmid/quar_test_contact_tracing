source("packages.R")
source("utils.R")
source("tracing_delays.R")
source("wolfel.R")
source("he.R")

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
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

run_analysis <- 
  function(n_sims          = 100,
           n_sec_cases     = 1000,
           n_ind           = 10000,
           seed            = 145,
           contact_info_delay,
           tracing_delay,
           asymp_parms){
    
    #browser()
    set.seed(seed)
    
    #Generate incubation periods to sample
    incubation_times <- make_incubation_times(
      n_travellers = n_ind,
      pathogen     = pathogen,
      asymp_parms = asymp_parms)
    
    inf <- tibble(sim=1:n_sims) %>%
      dplyr::mutate(prop_asy = rbeta(n = nrow(.),
                                     shape1 = asymp_parms$shape1,
                                     shape2 = asymp_parms$shape2)) 
    
    
    # Generate index cases' inc times
    ind_inc <- incubation_times %>% 
      filter(type=="symptomatic") %>% 
      sample_n(n_sims) %>% 
      mutate(sim=1:n_sims) %>% 
      left_join(inf) %>% 
      #add test delay (assume 2 days post onset)
      mutate(test_delay = 2) %>% 
      #sample contact info delay
      mutate(contact_info_delay=contact_info_delay %>% sample_n(1) %>% pull(t)) %>% 
      #sample tracing delay
      mutate(tracing_delay=tracing_delay %>% sample_n(1) %>% pull(t)) %>% 
      mutate(testing_t = onset+test_delay) %>% 
      mutate(traced_t = onset+test_delay+contact_info_delay+tracing_delay) %>% 
      mutate(sim=row_number()) 
        
    #Generate secondary cases
    sec_cases <- make_incubation_times(
      n_travellers = n_sec_cases,
      pathogen     = pathogen,
      asymp_parms = asymp_parms)
      
    ind_inc %<>% 
      nest(-c(sim,prop_asy,test_delay,contact_info_delay,tracing_delay,testing_t,traced_t)) %>% 
      mutate(sec_cases=map(.x=prop_asy,
                            incubation_times=sec_cases,
                            .f=make_sec_cases)) %>% 
      unnest(sec_cases)
    
    #exposure date relative to index cases exposure
    #filter out transmission prior to index case inf_start?
    ind_inc %<>% 
      mutate(exposed_t= si$r(n())) %>% 
      #obviously some better way to do this
      mutate(onset    = onset     + exposed_t,
             inf_start= inf_start + exposed_t,
             inf_end  = inf_end   + exposed_t,
             symp_end = symp_end  + exposed_t) %>% 
      #remove exposures post-test
      mutate(remove=ifelse(exposed_t>testing_t,TRUE,FALSE)) %>% 
      dplyr::filter(!remove) 
    
    #cross with scenarios
    incubation_times <- ind_inc %>% crossing(input) 
    
    source('kucirka_fitting.R',local=T)  
    
    #calc outcomes 
    incubation_times %<>% calc_outcomes(.,dat_gam)
    
    #when released
    incubation_times %<>% when_released()
    
    #stage of infection when released
    incubation_times %<>% stage_when_released()
    
    return(incubation_times)
    
  }

results <- run_analysis(contact_info_delay = getting_contact_info,
                        tracing_delay = tracing_delay,
                        asymp_parms = asymp_fraction)

results %>% 
  filter(stage_released=="Infectious") %>% 
  make_plots(.,input,faceting=~stringency,sum=F) 


