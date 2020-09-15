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

read(results_name) %>% 
    filter(yvar == "infectivity_averted",
           type == "all",
           #delay_scaling == 1,
           #stringency=="one",
           waning=="waning_none",
           index_test_delay == 2,
           quar_dur %in% c(0,14)
           ) %>%
    select(stringency, delay_scaling, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(iqr=paste0("(",iqr,")"),
           ui=paste0("(",ui,")")) %>% 
    select(delay_scaling,stringency,quar_dur,`50%`,iqr,ui) %>% 
    htmlTable()

# waning
read_results(results_name) %>% 
    filter(yvar == "infectivity_averted",
           type == "all",
           delay_scaling == 1,
           #stringency=="one",
           #waning=="waning_none",
           index_test_delay == 2,
           quar_dur %in% c(0,14)
    ) %>%
    select(stringency, waning, quar_dur, contains("%")) %>%
    mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
    unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
    unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
    mutate(iqr=paste0("(",iqr,")"),
           ui=paste0("(",ui,")")) %>% 
    select(stringency,waning,quar_dur,`50%`,iqr,ui) %>% 
    htmlTable()
