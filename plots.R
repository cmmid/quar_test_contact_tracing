
covid_pal2 <- c(All          = "black",
                Asymptomatic = "#3061ae",
                Symptomatic  = "#1DB954")

results_dat <- get(results_name) %>% 
  bind_rows() %>% 
  inner_join(input)

ribbon_plot <-
  function(x, 
           y_labels   = NULL, 
           colour_var = "stringency",
           faceting   = stringency ~ type){
    
    if (!any(grepl(pattern = "type", all.vars(faceting)))){
      x <- filter(x, type == "all")
    } 
    
    
    if (!is.null(y_labels)){
      x <- filter(x, yvar %in% names(y_labels))
      x <- mutate(x, yvar = factor(yvar, levels = names(y_labels), ordered = T))
      
    }
    
    # need to inser yvar
    
    f_lhs <- all.vars(lhs(faceting))
    f_rhs <- all.vars(rhs(faceting))
    
    faceting_new <-
      as.formula(
        paste(
          paste(f_lhs,
                collapse = " + "),
          paste(c("yvar", f_rhs),
                collapse = " + "),
          sep = " ~ "
        )
      )
    
    x %<>% mutate(stringency = capitalize(stringency),
                  type       = capitalize(type))
    
    colour_var_sym <- sym(colour_var)
    
    the_plot <- 
      ggplot(data = x, aes(x = second_test_delay,
                           y = M,
                           #color = !!colour_var,
                           fill  = !!colour_var_sym)) +
      facet_nested(nest_line = TRUE,
                   facets = faceting_new,
                   labeller = labeller(index_test_delay = index_test_labeller,
                                       delay_scaling    = delay_scaling_labeller,
                                       waning           = waning_labeller,
                                       stringency       = capitalize,
                                       type             = capitalize,
                                       yvar             = infectivity_labels)) +
      theme_minimal() +
      theme(legend.position = "bottom",
            panel.border = element_rect(fill=NA)) + 
      xlab("Time since exposure (days)") +
      ylab("Transmission potential") +
      geom_ribbon(aes(ymin = `2.5%`,
                      ymax = `97.5%`),
                  alpha = 0.2) +
      geom_ribbon(aes(ymin = `25%`,
                      ymax = `75%`),
                  alpha = 0.3) +
      geom_line(aes(y = `50%`,
                    color = !!colour_var_sym))
    
    
    if (colour_var == "stringency"){
      the_plot <- the_plot +
        scale_color_manual(name = "Number of tests required before release",
                           values = covid_pal) +
        scale_fill_manual(name = "Number of tests required before release",
                          values = covid_pal)
    } else{ 
      the_plot <- the_plot + 
        scale_color_manual(name = "Type of infection",
                           values = covid_pal2) +
        scale_fill_manual(name = "Type of infection",
                          values = covid_pal2)
    } 
    
    the_plot
  } 


ribbon_plot(results_dat,
            y_labels = infectivity_labels[c(3,2,1)],
            colour_var = "type",
            faceting = index_test_delay + type ~ stringency)

ggsave("results/ribbon_stringency.png", height=297/2, width=297,units="mm",dpi=400)
