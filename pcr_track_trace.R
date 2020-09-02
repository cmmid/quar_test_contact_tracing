source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("parameters.R")
source("kucirka_fitting.R")

results_name <- "test"

con <- file(paste0(results_name, ".log"))
sink(con, append=FALSE)
sink(con, append=TRUE, type="message")

input %<>% filter(stringency == "one",
                  index_test_delay  == 2,
                  second_test_delay == 0)

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

