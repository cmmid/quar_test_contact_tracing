results <- get(results_name)

results <- readRDS("results/sum_results.rds")


my_ribbon <- ribbon_plot(
  results_dat %>% filter(waning == "waning_none",
                         index_test_delay == 2,
                         delay_scaling == 1) %>% 
    #calculate time until release from exposure for each scenario
    mutate(time_since_exp=ifelse(stringency=="none",
                                 yes=quar_dur,
                                 no=quar_dur + results_delay * delay_scaling)),
  by_type = T,
  custom_facets = . ~ stringency ,
  y_labels = infectivity_labels["infectivity_averted"],
  ribbon = T,
  colour_var = "stringency")
main

ggsave("results/Figure_2_averted.png", width=210, height=210,units="mm",dpi=320)

delay_scaling <- ribbon_plot(
  results %>% 
    bind_rows() %>% 
    inner_join(input) %>% 
    filter(waning=="waning_none",
      index_test_delay==2,
      #delay_scaling == 1
    ) %>% 
    mutate(delay_scaling=fct_rev(as.factor(delay_scaling))) %>% View(),
  by_type = F,
  custom_facets =  ~  delay_scaling+stringency,
  y_labels = infectivity_labels[c("infectivity_averted")],
  ribbon = T,
  colour_var = "stringency")
delay_scaling

ggsave("results/delay_scaling.png", width=210, height=100,units="mm",dpi=320)
