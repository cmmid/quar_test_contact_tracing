
results_dat <- get(results_name) %>% 
  bind_rows() %>% 
  inner_join(input)


my_ribbon <- ribbon_plot(results_dat,
                         by_type = T,
                         y_labels = infectivity_labels[c(2,1)],
                         colour_var = "stringency")

ggsave("results/ribbon_stringency.png", height=210, width=297,units="mm",dpi=400)
