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
           index_test_delay    =  c(1,2,3),  # time to entering quarantine
           delay_scaling       =  c(1, 0.5)) %>%
  mutate(scenario=row_number()) 

con <- file("results.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

results <- 
  input %>%
  rowwise %>%
  group_split %>%
  map(~run_analysis(n_sims        = 1000,
                    n_ind_cases   = 1000,
                    n_sec_cases   = 1000,
                    seed          = 145,
                    input         = .x,
                    P_r           = P_r,
                    P_c           = P_c,
                    P_t           = P_t,
                    asymp_parms   = asymp_fraction))

results <- bind_rows(results, .id = "scenario") %>%
  mutate(scenario = parse_number(scenario)) %>%
  inner_join(input)

sink() 
sink(type="message")

source("figures_and_tables.R")
