results <- readRDS("results/sum_results_main_subset.rds")

results_dat <- results %>% 
  bind_rows() %>% 
  inner_join(input)

my_ribbon <- ribbon_plot(results_dat %>% filter(#waning=="waning_none",
                                                index_test_delay==2,
                                                delay_scaling == 1) %>% 
                           mutate(waning=fct_rev(waning)),
                         by_type = F,
                         #custom_facets = ~stringency,
                         y_labels = infectivity_labels["infectivity_avertable"],
                         ribbon=T,
                         colour_var = "stringency")
my_ribbon

ggsave("results/waning.png", width=210, height=100,units="mm",dpi=320)
