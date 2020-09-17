source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("parameters.R")
source("kucirka_fitting.R")


# input %<>% filter(
#   index_test_delay == 2,
#   delay_scaling    == 1,
#   waning           == "waning_canada_total",
#   #quar_dur         %in% c(0,5,7,10,14),
#   #stringency       == "none"
# )


input_split <-
  input %>% 
  rowwise %>%
  group_split

results_name <- "sum_results"

if (!dir.exists(here::here("results", results_name))){
  dir.create(here::here("results", results_name))
}

con <- file(here::here("results", results_name, "results.log"))
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
           return_full        = F
         )))

sink() 
sink(type="message")

saveRDS(get(results_name), 
        here::here("results", results_name, "results.rds"))
saveRDS(input,
        here::here("results", results_name, "input.rds"))

#source("figures.R")
source("plots.R")
source("tables.R")
source("sensitivities.R")

