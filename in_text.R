# in text numbers

## Reduction in transmission potential

results_index_test_delay_sensivity <- read_results("index_test_delay_sensivity")
    
results_index_test_delay_sensivity %>% 
    filter(yvar == "infectivity_pre",
           type == "all")

## longer quarantine reduces post-release 

results_sum_results_main_subset <- read_results("sum_results_main_subset")

results_sum_results_main_subset %>%
    filter(yvar == "infectivity_pre",
           type == "all",
           second_test_t %in% c(0, 7, 14)) # select on quar_dur

## changing TTI delays

