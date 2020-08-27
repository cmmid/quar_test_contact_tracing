results_dat <- get(results_name) %>% 
  bind_rows() 

results_dat %>% sample_n(100000) %>% 
  ggplot()+geom_jitter(aes(y=infectivity_averted,x=second_test_delay,colour=stringency),alpha=0.05)+
  facet_grid(stringency~type)+
scale_x_continuous(breaks=breaks_width(2))+
  scale_color_manual(name = "Number of negative tests required for release",
                     values = covid_pal)+
  theme(legend.position = "bottom")

ggsave("results/jitter.png",height=297,width=210,units="mm",dpi=600)
