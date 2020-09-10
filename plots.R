results <- readRDS("results/sum_results_main_subset.rds")

results_dat <- results %>% 
  bind_rows() %>% 
  inner_join(input)


my_ribbon <- ribbon_plot(
  results_dat %>% filter(waning=="waning_none",
                         index_test_delay==2,
                         delay_scaling == 1) %>% 
    #calculate time until release from exposure for each scenario
    mutate(time_since_exp=ifelse(stringency=="none",
                                 yes=quar_dur,
                                 no=quar_dur + results_delay * delay_scaling)),
  by_type = T,
  custom_facets = . ~ stringency ,
  y_labels = infectivity_labels[c(3)],
  ribbon = T,
  colour_var = "stringency")
my_ribbon

ggsave("results/no_ribbon.png", height=400, width=297,units="mm",dpi=600)
