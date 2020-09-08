
covid_pal2 <- c("#3061ae","#1DB954")
results_dat <- get(results_name) %>% 
  bind_rows() 

ribbon_plot <-
  function(x, 
           yvar = c("infectivity_averted" =
                      "Infectivity averted through quarantine and/or testing"), 
           colour_var = stringency,
           faceting = stringency ~ type,
           n = NULL, p = 1){
    
    if (!is.null(n)){x <- sample_n(x, n)} else {if (p < 1){x <- sample_frac(x, p)}}
    
    x$yvar <- x[[names(yvar)]]
    colour_var <- sym(colour_var)
    second_test_delays <- sort(unique(x$second_test_delay))
      
      the_plot <- 
        ggplot(data = x, aes(y = factor(second_test_delay, 
                                        levels = second_test_delays,
                                        ordered = T),
                             x = yvar,
                             color = !!colour_var,
                             fill  = !!colour_var)) +
        facet_grid(faceting) +
        theme_minimal() +
        theme(legend.position = "bottom",panel.border = element_rect(fill=NA)) + 
        xlab(yvar) +
        ylab("Time since exposure (days)") +
        coord_flip() +
        stat_summary(fun.min = function(x){
          set_names(quantile(x, probs = 0.25),
                    c("ymin"))},
          fun.max = function(x){
            set_names(quantile(x, probs = 0.75),
                      c("ymax"))},
          aes(group=!!colour_var), mult = 1, geom = 'ribbon',alpha=0.2,colour=NA)+
        stat_summary(fun.min = function(x){
          set_names(quantile(x, probs = 0.025),
                    c("ymin"))},
          fun.max = function(x){
            set_names(quantile(x, probs = 0.975),
                      c("ymax"))},
          aes(group=!!colour_var), mult = 1, geom = 'ribbon',alpha=0.2,colour=NA)+
        stat_summary(aes(group=!!colour_var), mult = 1, geom = 'smooth',alpha=0.3,
                     fun.data = function(x){
                       set_names(quantile(x, probs = c(0.025, 0.5, 0.975)),
                                 c("ymin", "y", "ymax"))})
      
         if (colour_var == "stringency"){
             the_plot <- the_plot + scale_color_manual(name = "Number of tests required before release",
                           values = covid_pal) +
                  scale_fill_manual(name = "Number of tests required before release",
                          values = covid_pal)
           } else{ 
             the_plot <- the_plot + scale_color_manual(name = "Type of infection",
                                          values = covid_pal2) +
                       scale_fill_manual(name = "Type of infection",
                                          values = covid_pal2)
        } 

    the_plot
  }

ribbon_plot(results_dat, n = 1350000, colour_var="stringency", faceting =  index_test_delay ~ stringency)

ggsave("results/ribbon_stringency.png",height=297,width=297,units="mm",dpi=400)
