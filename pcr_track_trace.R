source("packages.R")
source("utils.R")
source("tracing_delays.R")
#source("wolfel.R")
source("he.R")
source("kucirka_fitting.R")

waning_none <- function(x){
  waning_piecewise_linear(x, 1, 1, 7, 14)
}

waning_constant <- function(x){
  waning_piecewise_linear(x, 0.75, 0.75, 7, Inf)
}

waning_drop <- function(x){
  waning_piecewise_linear(x, 0.75, 0.25, 7, 14)
}


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
           delay_scaling       =  c(1, 0.5),
           waning              = "waning_drop") %>%
  mutate(scenario=row_number()) 

con <- file("results.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

input_split <- input %>%
  replace_na(replace = list(second_test_delay = 0)) %>%
  mutate(time_in_iso = as.numeric(screening) + first_test_delay + second_test_delay) %>%
  filter(time_in_iso %in% c(0, 3, 10, 14)) %>%
  rowwise %>%
  group_split

results_ <- 
  input_split %>%
  map(.x = .,
      ~run_analysis(n_sims        = 1000,
                    n_ind_cases   = 1000,
                    n_sec_cases   = 100,
                    input         = .x,
                    seed          = 145,
                    P_r           = P_r,
                    P_c           = P_c,
                    P_t           = P_t,
                    dat_gam       = dat_gam,
                    asymp_parms   = asymp_fraction))

sink() 
sink(type="message")

source("figures_and_tables.R")
