source("packages.R")
source("utils.R")
source("tracing_delays.R")
#source("wolfel.R")
source("he.R")
source("kucirka_fitting.R")

results_name <- "baseline_no_waning_w_sensitivities"

waning_none <- function(x){
  waning_points(x, X = 0, Y = 1)
}

waning_constant <- function(x){
  waning_points(x, X = 0, Y = 0.71)
}

# waning_drop <- function(x){
#   waning_piecewise_linear(x, 0.75, 0.25, 7, 14)
# }

# waning_linear <- function(x){
#   waning_piecewise_linear(x, ymax = 0.75, .16, 0, 8.3)
# }

waning_canada_community <- function(x){
  waning_points(x, X = c(0, 30), Y = c(1, 0.541), log = T)
}

waning_canada_total <- function(x){
  waning_points(x, X = c(0, 30), Y = c(1, 0.158), log = T)
}

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `none` = 
      crossing(screening         = FALSE,
               first_test_delay  = NA,
               second_test_delay = seq(0,14,by=2)), 
    `one` = 
      crossing(screening         = TRUE,
               first_test_delay  = NA,
               second_test_delay = seq(0,14,by=2)),
    `two` = 
      crossing(screening         = TRUE,
               first_test_delay  = 0,
               second_test_delay = seq(0,14,by=2))) %>% 
      bind_rows(.id = "stringency")) %>% 
  crossing(max_mip             = 14,
           post_symptom_window =  7,
           results_delay       =  1,
           index_test_delay    =  c(1, 2, 3),  # time to entering quarantine
           delay_scaling       =  c(1, 0.5),
           waning              = c("waning_none",
                                   "waning_constant",
                                   "waning_canada_total")) %>%
  filter(delay_scaling == 1, index_test_delay == 2, waning == "waning_none") %>%
  mutate(scenario=row_number()) 

con <- file(paste0(results_name, ".log"))
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

input_split <- input %>%
  #replace_na(replace = list(second_test_delay = 0)) %>%
  #mutate(time_in_iso = as.numeric(screening) + first_test_delay + second_test_delay) %>%
  #filter(time_in_iso %in% c(0, 3, 10, 14)) %>%
  rowwise %>%
  group_split

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f = ~run_analysis(
           n_sims        = 1000,
           n_ind_cases   = 1000,
           n_sec_cases   =  100,
           input         = .x,
           seed          =  145,
           P_r           = P_r,
           P_c           = P_c,
           P_t           = P_t,
           dat_gam       = dat_gam,
           asymp_parms   = asymp_fraction)))

sink() 
sink(type="message")

source("figures.R")
source("tables.R")
source("sensitivities.R")

