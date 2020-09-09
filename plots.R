
results_dat <- get(results_name) %>% 
  bind_rows() %>% 
  inner_join(input)


my_ribbon <- ribbon_plot(results_dat,
                         by_type = F,
                         y_labels = infectivity_labels[c(2,1)],
                         colour_var = "stringency")

ggsave("results/ribbon_stringency.png", height=297, width=297,units="mm",dpi=400)
