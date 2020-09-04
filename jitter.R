results_dat <- get(results_name) %>% 
  bind_rows() 

jitter_plot <- function(x, yvar = "infectivity_averted", 
                        faceting = stringency ~ type,
                        n = 1e5){

  x$yvar <- x[[yvar]]
  
  sample_n(x, n) %>%
    ggplot(data = ., aes(x = second_test_delay,
                         y = yvar,
                         fill = stringency)) +
    #geom_boxplot(,colour="black",outlier.alpha = 0)+
               position = position_jitter(w = 0.8, height = 0)) +
    facet_grid(faceting) +
    scale_x_continuous(breaks=breaks_width(2))+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme(legend.position = "bottom") + 
    ylab(yvar)
  
}

jitter_plot(results_dat, n = 5e4, yvar = "infectivity_post")

ggsave("results/jitter.png",height=297,width=210,units="mm",dpi=600)
