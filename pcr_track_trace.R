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

input %<>% filter(index_test_delay == 2)

input_split <-
  input %>% 
  rowwise %>%
  group_split

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f = ~run_analysis(
           n_sims             = 100,
           n_ind_cases        = 1000,
           n_sec_cases        =  100,
           input              = .x,
           seed               =  145,
           P_r                = P_r,
           P_c                = P_c,
           P_t                = P_t,
           dat_gam            = dat_gam,
           asymp_parms        = asymp_fraction,
           return_full        = FALSE,
           y_labels           = infectivity_labels
           )))

sink() 
sink(type="message")

#source("figures.R")
source("tables.R")
source("sensitivities.R")

