## functions used for ploting

# colours
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
covid_pal2 <- set_names(covid_pal, c("All", "Asymptomatic", "Symptomatic"))
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)

infectivity_labels <-
  c("infectivity_post" =
      "Transmission potential of secondary cases \nafter release",
    "infectivity_averted" = 
      "Transmission potential of secondary cases \naverted as a result of quarantine and testing",
     "infectivity_quar" = 
       "Transmission potential in community due to imperfect quarantine adherence",
    "infectivity_pre" =
      "Transmission potential of secondary cases \nprior to being traced",
    "infectivity_total" = 
      "Transmission potential of secondary cases \nin community compared to no quarantine or testing"
  )



test_labeller <- function(x){
  mutate(x,
         stringency = factor(stringency,
                             levels = c("none",
                                        "one",
                                        "two"),
                             labels = c("None",
                                        "One",
                                        "Two"),
                             ordered = T))
}

type_labels <- c("asymptomatic" =
                   "Asymptomatic",
                 "symptomatic" =
                   "Pre-symptomatic")

index_test_labeller <- function(x, newline = FALSE){
  paste0("Delay between index case's\nonset and having a test:",
         ifelse(newline, "\n", " "),
         x,
         ifelse(x == 1, " day", " days"))
}

delay_scaling_labeller <- function(x, newline = FALSE){
  paste0("Contact tracing delay",
         ifelse(newline, "\n", " "),
         "scaling factor: ",
         x)
}


waning_labeller <- function(x){
  paste("Adherence to quarantine guidance:\n",
        dplyr::case_when(
          x == "waning_canada_total"     ~ "Exponential decay",
          x == "waning_constant"         ~ "Constant",
          x == "waning_none"             ~ "No waning",
          x == "waning_canada_community" ~ "Exponential decay (community only)",
          TRUE ~ "Unknown"))
}

percentage <- function(x, ...){
  
  ans <- scales::percent(x, ...)
  ans[x > 1] <- ""
  
  ans
  
}

pretty_percentage <- function(x){
  ans <- pretty(x)
  ans[ans <= 1]
}


make_release_figure <- function(x_summaries,
                                input,
                                xlab = "Days in quarantine",
                                ylab = "",
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                hline = 0,
                                faceting = NULL,
                                percent = FALSE){
  
  x_summaries %<>% test_labeller # should this be in the facet call?
  
  # how to do presymptomatic
  
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x_summaries %<>% mutate(type = factor(type,
                                          levels = c("asymptomatic",
                                                     "symptomatic"),
                                          labels = c("Asymptomatic",
                                                     "Presymptomatic")))
  }
  
  
  figure <-  
    ggplot(data=x_summaries, aes(x = time_since_exp, 
                                 y = `50%`, 
                                 color = stringency)) +
    geom_hline(aes(yintercept=1), linetype=hline)+
    geom_linerange(aes(ymin  = `2.5%`, 
                       ymax  = `97.5%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.3,
                   size      = 3) +
    geom_linerange(aes(ymin  = `25%`,
                       ymax  = `75%`,
                       group = stringency),
                   position  = position_dodge2(width = 1),
                   alpha     = 0.5,
                   size      = 3) +
    geom_point(pch           = "-",
               size          = 12,
               position      = position_dodge2(width = 1),
               aes(y         = `50%`,
                   group     = stringency)) +
    scale_x_continuous(breaks = breaks_width(2))+
    scale_color_manual(name = "Number of negative tests required for release",
                       values = covid_pal)+
    theme_minimal()+
    theme(axis.ticks = element_line(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.position = "bottom",
          strip.placement = "outside") +
    ggplot2::labs(x = xlab,
                  y = ylab) +
    xlab("Days in quarantine\n(including 1 day delay on testing results)")
  
  figure <- figure + 
    facet_nested(nest_line = TRUE,
                 facets = faceting,
                 labeller = labeller(index_test_delay = index_test_labeller,
                                     delay_scaling    = delay_scaling_labeller,
                                     waning           = waning_labeller))
  
  
  return(figure)
  
}

plot_data <- function(input, 
                      x_summaries,
                      main_scenarios = NULL){
  
  dat <- x_summaries  %>%
    inner_join(input) %>% # should really carry these through when summarising
    mutate(time_since_exp_ = 
             ifelse(is.na(time_since_exp),
                    0,
                    time_since_exp),
           time_in_iso = 
             first_test_delay + 
             time_since_exp_+
             screening)
  
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% dplyr::select(-one_of("released_test")) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    #tidyr::nest(data = -c(first_test_delay, time_since_exp)) %>%
    dplyr::mutate(delays = paste(first_test_delay, "&",
                                 first_test_delay + time_since_exp)) %>%
    #tidyr::unnest(data) %>%
    dplyr::mutate(time_in_iso = factor(time_in_iso, 
                                       levels = sort(unique(.$time_in_iso)),
                                       ordered = T)) %>%
    dplyr::filter(M!=0) %>%  # if the mean is zero, this group is empty
    return
}

make_arrivals_table <- function(x, table_vars = c("country")){
  x %>%
    mutate(sim = factor(sim, levels = 1:n_arrival_sims)) %>%
    group_by_at(.vars = vars(sim, one_of(table_vars))) %>%
    count(.drop = F) %>%
    group_by_at(.vars = vars(-n, -sim)) %>%
    nest() %>%
    mutate(Q = map(data, ~quantile(.x$n, probs = probs))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    group_by_at(.vars = vars(one_of(table_vars))) %>%
    transmute(value = sprintf("%0.0f (95%%: %0.0f, %0.0f)", `50%`, `2.5%`, `97.5%`)) %>%
    mutate_at(.vars = vars(-c(country, value)), 
              .funs = str_to_title) %>%
    spread(country, value)
}


save_plot <- function(plot   = ggplot2::last_plot(),
                      prefix = stringi::stri_rand_strings(1, length = 8),
                      base   = NULL, # to help identify which figure in the paper
                      device = NULL, # additional devices, e.g. "png" 
                      width  = 210, 
                      height = 210,
                      dpi    = 300,
                      units  = "mm"){
  
  file <- paste0("results/",
                 paste(prefix, base, sep = "_"),
                 ".pdf")
  
  pdf(file = file,
      width  = measurements::conv_unit(x = width,  units, "inch"),
      height = measurements::conv_unit(x = height, units, "inch"),
      useDingbats = FALSE)
  
  print(plot)
  
  dev.off()
  
  if (length(device) > 0L){
    #img <- pdftools::pdf_render_page(file, dpi = dpi)
    
    img <- magick::image_read_pdf(file, density = dpi)
    
    purrr::map(.x = device,
               .f = 
                 ~magick::image_convert(image = img, format = .x) %>%
                 magick::image_write(., path = sub(pattern = "pdf",
                                                   replacement = .x,
                                                   x = file)))
  }
}



make_days_plots <- 
  function(x, 
           main_scenarios   = NULL,
           plot             = TRUE,
           log_scale        = FALSE,
           text_size        = 2.5,
           xlab             = "Days since exposure\n(including 1 day delay on testing results)",
           sum              = FALSE,
           y_labels         = NULL, # pass in y_vars as a named list
           faceting         = NULL,
           dir              = stringi::stri_rand_strings(1, 8),
           base             = stringi::stri_rand_strings(1, 8)){
    
    if (!dir.exists(paste0("results/",dir))){
      dir.create(paste0("results/",dir))
    }
    
    all_grouping_vars <- all.vars(faceting)
    
    if (sum){
      y_labels <- sub("^Average", "Total", y_labels)
    } 
    
    x_days_summaries <-
      as.list(names(y_labels)) %>%
      lapply(X = ., FUN = function(y){
        map_df(.x = x,
               ~make_released_time_quantiles(.x,
                                             y_var = y, sum = sum,
                                             vars = all_grouping_vars))})
    ## end summaries
    
    fig_data <- x_days_summaries %>% 
      map(~plot_data(input = input, # should we still pass it in? recover?
                     x_summaries = 
                       .x,
                     main_scenarios))
    if (plot){
      
      figs <- map2(
        .x = fig_data,
        .y = y_labels,
        .f = ~make_release_figure(
          x         = .x,
          #input     = input,
          xlab      = xlab,
          text_size = text_size,
          ylab      = .y,
          faceting  = faceting,
          percent   = TRUE) )
      
      
      if (length(figs) > 1){
        legend <- cowplot::get_legend(figs[[1]])
        figs   <- lapply(X = figs, function(x){x + theme(legend.position = "none")})
        
        fig    <- cowplot::plot_grid(plotlist = figs, nrow = 1,
                                     labels = LETTERS[1:length(figs)])
        
        fig    <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .2))
        
      } else {
        fig <- figs[[1]]
      }
      
      
      list("png", "pdf") %>%
        map(~ggsave(filename = paste0("results/",dir,"/days_plots_",base,".",.x),
                    plot=fig,
                    width  = 60*nrow(distinct(fig_data[[1]][,get.vars(rhs(faceting))]))*
                      length(fig_data), 
                    height = 120*nrow(distinct(fig_data[[1]][,get.vars(lhs(faceting))])), 
                    dpi = 600,
                    units = "mm",
                    device = ifelse(.x == "pdf",
                                    cairo_pdf,
                                    "png")))
    }
    
    return(
      set_names(fig_data, sub(pattern = "infectivity_", 
                              replacement = "", x = names(y_labels)))
      
    )
    
  }


summarise_results <- function(x, reduction = TRUE){
  if (!is.logical(reduction)){
    stop("reduction must be logical")
  }
  x <- mutate_at(x, 
                 .vars = vars(contains("%")),
                 .funs = function(x){
                   reduction*(1 - x) +
                     (1 - reduction)*x})
  if (reduction){
    percentages <- grep(x = names(x), pattern = "%")
    x[,rev(percentages)] <- x[,percentages]  
  }
  
  x
  
}

show_results <- function(x, reduction = TRUE){
  dplyr::select(x, delays,
                one_of(all.vars(faceting)),
                screening, time_in_iso,
                contains("%")) %>%
    group_by(stringency, index_test_delay) %>%
    group_split %>%
    map(summarise_results, reduction = reduction) 
}



ribbon_plot <-
  function(x, 
           y_labels   = NULL, 
           colour_var = "stringency",
           by_type    = FALSE,
           custom_facets = NULL,
           ribbon  =TRUE
  ){

    
    if (is.null(custom_facets)){
      f_lhs <- c("waning", "index_test_delay", "delay_scaling")
    f_rhs <- c("yvar", "stringency")
    
    } else {
      f_lhs <- all.vars(lhs(custom_facets))
      f_rhs <- all.vars(rhs(custom_facets))
    }
    
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
    
    xlims <- range(x$time_since_exp)
    
    colour_var_sym <- sym(colour_var)
    
    the_plot <- 
      ggplot(data = x, aes(x = time_since_exp,
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
      ylab("Transmission potential") 
    
    if(ribbon == TRUE){
     the_plot <-  the_plot + 
       geom_ribbon(aes(ymin = `2.5%`,
                      ymax = `97.5%`),
                  alpha = 0.2) +
      geom_ribbon(aes(ymin = `25%`,
                      ymax = `75%`),
                  alpha = 0.3)
    } else {
       the_plot <- the_plot +
         geom_line(aes(y=`2.5%`,colour=!!colour_var_sym),linetype="dashed")+
         geom_line(aes(y=`97.5%`,colour=!!colour_var_sym),linetype="dashed")
     }
    
    the_plot <- the_plot +
      geom_line(aes(y = `50%`,
                    color = !!colour_var_sym)) +
      scale_x_continuous(minor_breaks = seq(xlims[1], xlims[2], by = 1),
                         breaks       = seq(xlims[1], xlims[2], by = 2))
    
    
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


