source('run_analysis_func.R')

test <- run_analysis(
  n_arrival_sims  = 100,
  countries       = c("EU", "USA"),
  trav_vol_manual = 10000,
  fixed           = TRUE,
  n_travellers    = 1e4,
  trav_vol_p      = 1,
  seed            = 145)

test_fig <-
  make_plots(arrival_released_times = test, 
             input = input,
             text_size = 2.5,
             log_scale = FALSE,
             main_scenarios = main_scenarios,
             faceting = country ~ stringency,
             pre_board_screening = FALSE)

save_plot(test_fig, width = 210, height = 210, units = "mm", device = NULL)

test_fig_type <-
  make_plots(arrival_released_times = test, 
             input = input,
             text_size = 2.5,
             log_scale = FALSE,
             main_scenarios = main_scenarios,
             faceting = country + type ~ stringency,
             pre_board_screening = FALSE)

save_plot(test_fig_type, base = "type", 
          height = 420, width = 210,
          device = c("pdf", "png"))

test_real <- run_analysis(
  n_arrival_sims  = 1000,
  countries       = c("EU", "USA"),
  n_travellers    = 1e4,
  trav_vol_p      = 7/30,
  fixed           = FALSE,
  seed            = 145)
