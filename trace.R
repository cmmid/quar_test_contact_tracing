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
           n_sec_cases     = 10,
           seed            = 145,
           contact_info_delay,
           tracing_delay,
           asymp_parms){
  
    #browser()
    set.seed(seed)
    

    #gen index case
    index <- make_incubation_times(n_travellers = n_sims,
                          pathogen = pathogen,
                          asymp_parms = asymp_parms) %>% 
      filter(type=="symptomatic") %>% 
      rename_all(paste0,"_index") %>% 
      #add test delay (assume 2 days post onset)
      mutate(test_delay = 2) %>% 
      #sample contact info delay
      mutate(contact_info_delay=contact_info_delay %>% sample_n(1) %>% pull(t)) %>% 
      #sample tracing delay
      mutate(tracing_delay=tracing_delay %>% sample_n(1) %>% pull(t)) %>% 
      mutate(testing_t = onset_index+test_delay) %>% 
      mutate(traced_t = onset_index+test_delay+contact_info_delay+tracing_delay) %>% 
      rename("sim"=idx_index) %>% 
      group_by(sim) %>% 
      nest(.key="index_case_data") 
      

    
   #do secondary cases need to be produced by negative binomial?
   sec_cases <- make_incubation_times(n_travellers = n_sec_cases,
                         pathogen      = pathogen,
                         asymp_parms = asymp_parms) %>% 
     sample_n(n_sec_cases) %>% 
     nest()
  
   #join secondary cases to index
    transmission <- index %>% 
      mutate(sec_cases=sec_cases$data) %>% 
      unnest(index_case_data) %>% 
      unnest(sec_cases)
  
    
    #exposure date relative to index cases exposure
    #filter out transmission prior to index case inf_start?
    transmission %<>% 
      mutate(exposed_t= si$r(n())) %>% 
      #obviously some better way to do this
      mutate(onset    = onset     + exposed_t,
             inf_start= inf_start + exposed_t,
             inf_end  = inf_end   + exposed_t,
             symp_end = symp_end  + exposed_t)
    
    
    #remove exposures post-test
    transmission %<>% 
      mutate(remove=ifelse(exposed_t>testing_t,TRUE,FALSE)) %>% 
      dplyr::filter(!remove) 
    
    
    #cross with scenarios
    incubation_times <- transmission %>% crossing(input) 
    
    source('kucirka_fitting.R',local=T)  
    
    #calc outcomes 
    incubation_times %<>% calc_outcomes(.,dat_gam)
    
    #when released
    incubation_times %<>% when_released()
    
    incubation_times %<>% stage_when_released()
    
    return(incubation_times)
    
  }

results <- run_analysis(contact_info_delay = getting_contact_info,
                        tracing_delay = tracing_delay,
                        asymp_parms = asymp_fraction)

results %>% 
  filter(stage_released=="Infectious") %>% 
  inner_join(input) %>% 
  make_plots(.,input,faceting=~stringency) 


