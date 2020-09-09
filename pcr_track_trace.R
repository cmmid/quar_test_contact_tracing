source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("parameters.R")
source("kucirka_fitting.R")


input %<>% filter(
  index_test_delay == 2,
    delay_scaling  == 1,
    waning         == "waning_none",
    stringency     == "one")

nrow(input)

input_split <-
  input %>% 
  rowwise %>%
  group_split

results_name <- "test"
con <- file(paste0(results_name, "2.log"))
sink(con, append=FALSE)
sink(con, append=TRUE, type="message")

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f = ~run_analysis(
           n_sims             = 1000,
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

saveRDS(get(results_name),"results/sum_results_subset.rds")

#source("figures.R")
source("plots.R")
source("tables.R")
source("sensitivities.R")

