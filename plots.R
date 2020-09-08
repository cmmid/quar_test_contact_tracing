
covid_pal2 <- set_names(covid_pal, c("All", "Asymptomatic", "Symptomatic"))

results_dat <- get(results_name) %>% 
  bind_rows() %>% 
  inner_join(input)

ribbon_plot <-
  function(x, 
           y_labels   = NULL, 
           colour_var = "stringency",
           by_type    = FALSE
           ){
    
    f_lhs <- c("waning", "index_test_delay", "delay_scaling")
    f_rhs <- c("yvar", "stringency")
    
    if (!by_type){
      x <- filter(x, type == "all")
    } else {
     f_lhs <- c(f_lhs, "type") 
    }
    

    if (!is.null(y_labels)){
      x <- filter(x, yvar %in% names(y_labels))
      x <- mutate(x, yvar = factor(yvar, levels = names(y_labels), ordered = T))
    }
    
    # here we want to drop anything that's only got one value
    
    drop_unused_terms <- function(y,x){
      map_chr(.x = y, .f = function(v){
        if (length(unique(x[[v]])) > 1L){
          v
        } else {NA_character_}
      }) %>% na.omit %>% c
    }
    
    f_lhs <- drop_unused_terms(f_lhs, x)
    f_rhs <- drop_unused_terms(f_rhs, x)
    
    
    faceting_new <-
      as.formula(
        paste(
          paste(f_lhs,
                collapse = " + "),
          paste(f_rhs,
                collapse = " + "),
          sep = " ~ "
        )
      )
    
    x %<>% mutate(stringency = capitalize(stringency),
                  type       = capitalize(type))
    
    xlims <- range(x$second_test_delay)
    
    colour_var_sym <- sym(colour_var)
    
    the_plot <- 
      ggplot(data = x, aes(x = second_test_delay,
                           y = M,
                           #color = !!colour_var,
                           fill  = !!colour_var_sym)) +
      facet_nested(nest_line = TRUE,
                   drop      = TRUE,
                   facets    = faceting_new,
                   labeller  = labeller(index_test_delay = index_test_labeller,
                                        delay_scaling    = 
                                          function(x){delay_scaling_labeller(x,TRUE)},
                                        waning           = waning_labeller,
                                        stringency       = capitalize,
                                        type             = capitalize,
                                        yvar             = infectivity_labels,
                                        .multi_line      = TRUE)) +
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
                    color = !!colour_var_sym)) +
      scale_x_continuous(minor_breaks = seq(xlims[1], xlims[2], by = 1),
                         breaks       = seq(xlims[1], xlims[2], by = 7))
    
    
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


my_ribbon <- ribbon_plot(results_dat,
                         by_type = F,
                         y_labels = infectivity_labels[c(2,1)],
                         colour_var = "stringency")

ggsave("results/ribbon_stringency.png", height=297, width=297,units="mm",dpi=400)
