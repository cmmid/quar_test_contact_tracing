source("packages.R")
source("utils.R")
source("tracing_delays.R")
#source("wolfel.R")
source("he.R")

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `low` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = 0,
               second_test_delay = NA), 
    `moderate` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = c(3,5,7,9),
               second_test_delay = NA),
    `high` = 
      crossing(screening         = TRUE,
               first_test_delay  = c(0:3),
               second_test_delay = c(2,4,6)),
    `maximum` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = 14,
               second_test_delay = NA)) %>%
      bind_rows(.id = "stringency")) %>% 
  crossing(max_mip             = 14,
           post_symptom_window =  7,
           results_delay       =  1,
           index_test_delay    =  c(1,2,3)) %>% # time to entering quarantine
  mutate(scenario=row_number()) 

results <- run_analysis(n_sims             = 1000,
                        n_ind_cases        = 10000,
                        n_sec_cases        = 1000,
                        contact_info_delay = getting_contact_info,
                        index_result_delay = index_result_delay,
                        tracing_delay      = tracing_delay,
                        asymp_parms        = asymp_fraction)

source("figures_and_tables.R")
