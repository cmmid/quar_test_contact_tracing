#produce plots and figures for text
source("run_analysis_func.R")

run_per_10000_analysis(text_size = 2.5)

gc()
real_trav_vol <- run_real_trav_vol_analysis()
real_trav_vol$table_1
real_trav_vol$main_text_dat

write_csv(real_trav_vol$n_data, "results/trav_vol_n_results.csv")
write_csv(real_trav_vol$pd_data,"results/trav_vol_pd_results.csv")

gc()
# i feel like these should be specifications of relevant scenario properties
# scenario 7 might change with different parameter values in setup
RR_7 <- run_rr_analysis(baseline_scenario = "7",log_scale = FALSE,trav_vol_manual = 100000,n_arrival_sims = 500)
write.csv(RR_7$n_fig_data,"results/low_baseline_n_RR_results.csv")
write.csv(RR_7$pd_fig_data,"results/low_baseline_person_days_RR_results.csv")

gc()
RR_87 <- run_rr_analysis(baseline_scenario = "87",log_scale = TRUE,trav_vol_manual = 100000,n_arrival_sims = 500)
write.csv(RR_87$n_fig_data,"results/max_baseline_n_RR_results.csv")
write.csv(RR_87$pd_fig_data,"results/max_baseline_person_days_RR_results.csv")

# come on! what's this!?
real_trav_vol$n_data %>% filter(first_test_delay==7,!post_flight_screening,pre_board_screening=="None") 
real_trav_vol$n_data %>% filter(first_test_delay==14,!post_flight_screening,pre_board_screening=="None") 

RR_7$n_fig_data %>% filter(first_test_delay==7,!post_flight_screening,pre_board_screening=="None") 
RR_7$n_fig_data %>% filter(first_test_delay==14,!post_flight_screening,pre_board_screening=="None") 

RR_7$pd_fig_data %>% filter(first_test_delay==7,!post_flight_screening,pre_board_screening=="None") %>% mutate_at(vars(`2.5%`:`97.5%`),function(x){(1-x)*100})
RR_7$pd_fig_data %>% filter(first_test_delay==14,!post_flight_screening,pre_board_screening=="None") %>% mutate_at(vars(`2.5%`:`97.5%`),function(x){(1-x)*100})

RR_87$n_fig_data %>% filter(first_test_delay==3,second_test_delay==4,pre_board_screening=="None") 

RR_7$n_fig_data %>% filter(first_test_delay==7,post_flight_screening,pre_board_screening=="None") 

source("kucirka_plot.R")