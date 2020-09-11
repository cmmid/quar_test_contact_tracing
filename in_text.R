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
           quar_dur %in% c(0, 7, 14)) # select on quar_dur rather than time_since_exp

## changing TTI delays

results_delay_scaling_sensivity <- read_results("delay_scaling_sensivity")

results_delay_scaling_sensivity %>%
    filter(yvar == "infectivity_averted",
           type == "all",
           index_test_delay == 2,
           quar_dur %in% c(0, 7, 14)) %>%
    select(stringency, delay_scaling, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1))
