source("packages.R")
source("utils.R")
source("tracing_delays.R")
source("wolfel.R")
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
               first_test_delay  = c(3,5,7),
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
  mutate(scenario=row_number()) %>%
  mutate(max_mqp             = 14,
         post_symptom_window =  7,
         results_delay       =  2,
         test_delay          =  2
  )


results <- run_analysis(contact_info_delay = getting_contact_info,
                        tracing_delay      = tracing_delay,
                        asymp_parms        = asymp_fraction)

results %>% 
  filter(stage_released=="Infectious") %>% 
  make_plots(.,input, 
             faceting = ~stringency,
             y_var = "days_prior_inf",
             sum = F) 


