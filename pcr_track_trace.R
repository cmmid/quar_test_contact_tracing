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
  crossing(max_mqp             = 14,
           post_symptom_window =  7,
           results_delay       =  1,
           index_test_delay    =  c(0,2,4)) %>% # time to entering quarantine
  mutate(scenario=row_number()) 

results <- run_analysis(n_sims=1000,
                        n_sec_cases = 1000,
                        contact_info_delay = getting_contact_info,
                        index_result_delay = index_result_delay,
                        tracing_delay      = tracing_delay,
                        asymp_parms        = asymp_fraction)

results %>% make_plots(.,input, 
                            faceting = index_test_delay ~ stringency,
                            y_var = "days_released_inf",
                            sum = F)

results_df <- results %>% make_days_plots(.,input, 
                            faceting = index_test_delay ~ stringency,
                            y_vars = c("days_prior_inf","days_released_inf"),
                            sum = F)

baseline_low <- data.frame(
  screening = FALSE,
  first_test_delay      = 0,
  second_test_delay     = NA,
  stringency            = "low"
)


low <- run_rr_analysis(results,
                      main_scenarios, 
                      baseline_scenario = baseline_low,
                      faceting = index_test_delay ~ stringency,
                      log_scale=F)

low %>% filter(time_in_iso>8,index_test_delay==2) %>% mutate_at(vars(`2.5%`:`97.5%`),function(x) 1- x)

baseline_max <- data.frame(
  screening = FALSE,
  first_test_delay      = 14,
  second_test_delay     = NA,
  stringency            = "maximum"
)
max <- run_rr_analysis(results,
                       main_scenarios, 
                       baseline_scenario = baseline_max,
                       faceting = index_test_delay ~ stringency,
                       log_scale=TRUE)




