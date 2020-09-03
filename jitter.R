results_dat <- get(results_name) %>% 
  bind_rows() 

jitter_plot <- function(x, yvar = "infectivity_averted", 
                        geom = c("point"),
                        faceting = stringency ~ type,
                        n = 1e5, p = 0.1){
  
  if (!is.null(n)){
    x <- sample_n(x, n)
  } else {
    x <- sample_frac(x, p)
  }
  
  
  x$yvar <- x[[yvar]]
  second_test_delays <- sort(unique(x$second_test_delay))
  
  the_plot <- 
    ggplot(data = x, aes(x = factor(second_test_delay, 
                                    levels = second_test_delays,
                                    ordered = T),
                         y = yvar,
                         color = stringency)) +
    facet_grid(faceting) +
    #scale_x_continuous(breaks=breaks_width(2))+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme(legend.position = "bottom") + 
    ylab(yvar)
  
  if ("point" %in% geom){
    the_plot <- the_plot + 
      geom_point(alpha = 0.05, pch = 16,
                 position = position_jitter(w = 0.8, height = 0))
  }
  
  if ("violin" %in% geom){
    the_plot <- the_plot + 
      geom_half_violin(aes(fill = stringency), 
                       side = "r",
                       draw_quantiles = probs,
                       alpha = 0.5) +
      geom_half_boxplot(side = "l") +
      scale_fill_manual(name = "Number of negative tests required for release",
                        values = covid_pal)
  }
  
  if ("hist" %in% geom){
    
    nudge_x_value <- -0.15
    
    the_plot <- 
      ggplot(data = x, aes(y = factor(second_test_delay, 
                                      levels = second_test_delays,
                                      ordered = T),
                           x = yvar,
                           color = stringency,
                           fill  = stringency)) +
      facet_grid(faceting) +
      geom_density_ridges(stat = "binline", 
                          color = NA,
                          panel_scaling = FALSE,
                          #aes(y = ..count..),
                          breaks = seq(0,1,by = 0.05), 
                          pad = FALSE,
                          scale = 0.9 + nudge_x_value) +
      #scale_x_continuous(breaks=breaks_width(2))+
      scale_color_manual(name = "Number of negative tests required for release",
                         values = covid_pal) +
      scale_fill_manual(name = "Number of negative tests required for release",
                        values = covid_pal) + 
      theme_bw() +
      theme(legend.position = "bottom") + 
      xlab(yvar) +
      ylab("Duration of quarantine (days)") +
      coord_flip() +
      stat_summary(position = position_nudge(x = 0, y = nudge_x_value),
                   geom = "linerange",
                   size = 2,
                   alpha = 0.25,
                   fun.data = function(x){
                     set_names(quantile(x, probs = c(0.025, 0.5, 0.975)),
                               c("ymin", "y", "ymax"))}
      ) +
      
      stat_summary(position = position_nudge(x = 0, y = nudge_x_value),
                   geom = "linerange",
                   size = 2,
                   alpha = 0.25,
                   fun.data = function(x){
                     set_names(quantile(x, probs = c(0.25, 0.5, 0.75)),
                               c("ymin", "y", "ymax"))}
      ) +
      
      stat_summary(position = position_nudge(x = 0, y = nudge_x_value),
                   geom = "pointrange",
                   fill = "white",
                   solid = TRUE,
                   size = 0.5,
                   shape = 21,
                   stroke = 0.5,
                   fun.data = function(x){c(ymin = NA,
                                            y    = mean(x),
                                            ymax = NA)}) +
      stat_summary(position = position_nudge(x = 0, y = nudge_x_value),
                   geom = "pointrange",
                   size = 2,
                   pch = "-",
                   fun.data = function(x){c(ymin = NA,
                                            y    = median(x),
                                            ymax = NA)})  
    
  } 
  
  the_plot
}

jitter_plot(results_dat, n = 5e4, faceting = stringency ~ .,
            yvar = "infectivity_averted", geom = c("hist"))

ggsave("results/jitter.png",height=297,width=210,units="mm",dpi=600)
