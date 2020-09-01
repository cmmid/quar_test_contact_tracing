## functions used for ploting

# colours
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)


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
        dplyr::case_when(x == "waning_canada_total" ~ "Exponential decay",
                         x == "waning_constant"     ~ "Constant",
                         x == "waning_none"         ~ "No waning",
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
  #browser()
  x_summaries %<>% test_labeller
  
  # how to do presymptomatic
  
  facet_vars <- all.vars(faceting)
  
  if ("type" %in% facet_vars){
    x_summaries %<>% mutate(type = factor(type,
                                          levels = c("asymptomatic",
                                                     "symptomatic"),
                                          labels = c("Asymptomatic",
                                                     "Presymptomatic")))
  }
  
  require(formula.tools)
  
  
  # dy <- dplyr::select(x_summaries, !!!facet_vars,
  #                     `97.5%`, `75%`) %>%
  #   gather(key, value, -c(!!!facet_vars)) %>%
  #   group_by_at(.vars = vars(-value)) %>%
  #   filter(value == max(value)) %>%
  #   distinct %>%
  #   spread(key, value) %>%
  #   mutate(stringency = fct_collapse(stringency, 
  #                                    "Two" = "two",
  #                                    other_level = "Other")) %>%
  #   group_by_at(.vars = vars(-`75%`, -`97.5%`)) %>%
  #   summarise_all(.funs = ~max(.)) %>%
  #   tidyr::pivot_wider(names_from = "stringency", 
  #                      names_glue = "{stringency}_{.value}",
  #                      values_from = c(`75%`, `97.5%`))
  
  
  
  # needs_expanding <- dy %>% 
  #   filter(`High_97.5%` > 0.75*`Other_97.5%`) %>%
  #   {nrow(.) > 0L}
  # 
  # dy %<>% mutate(ypos = 0.05 * pmax(`High_97.5%`, `Other_97.5%`)) %>%
  #   dplyr::select(one_of(all.vars(faceting)), ypos)
  
  
  
  #x_summaries %<>% left_join(dy)
  
  figure <-  
    ggplot(data=x_summaries, aes(x = second_test_delay, 
                                 y = `50%`, 
                                 color = stringency)) +
    geom_hline(aes(yintercept=1),linetype=hline)+
    geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`,
                       group = stringency),
                   position = position_dodge2(width = 1),
                   alpha = 0.3,
                   size = 3) +
    geom_linerange(aes(ymin = `25%`, ymax = `75%`,
                       group = stringency),
                   position = position_dodge2(width = 1),
                   alpha = 0.5,
                   size = 3) +
    geom_point(pch = "-", size = 12,
               position = position_dodge2(width = 1),
               aes(y = `50%`,
                   group = stringency)
    ) +
    # geom_text(data=filter(x_summaries, stringency=="High"),
    #           aes(x     = time_in_iso,
    #               y     = `97.5%` + ypos,#ifelse(percent, 1.05, `97.5%` + ypos),
    #               label = delays),
    #           angle     = text_angle,
    #           hjust     = h_just,
    #           vjust     = 0,
    #           size      = text_size,
    #           position  = position_dodge2(width = 0.75),
    #           check_overlap = TRUE,
    #           show.legend = F)+
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
  
  # figure <- figure+ 
  #     facet_nested(
  #   nest_line = T,
  #   facets = faceting,
  #   labeller = labeller(index_test_delay = index_test_labeller,
  #                       delay_scaling    = delay_scaling_labeller,
  #                       waning           = waning_labeller),
  #   scales = "free_x", space = "free")
  # }
  
  # check if the top of the y axis needs adjustment
  
  
  
  # end check top y
  
  # if (log_scale) {
  #   mult <- c(0.1, ifelse(needs_expanding, 0.5, 0.1))
  #   figure <- figure +
  #     scale_y_log10(labels=format_format(scientific=FALSE),
  #                   expand = expansion(mult = mult)) +
  #     theme(panel.grid.minor.y = element_blank()) +
  #     annotation_logticks(sides = "l")
  #   return(figure)
  #   
  # } 
  # 
  # if (percent){
  #   mult <- c(0.01, ifelse(needs_expanding, 0.25, 0.1))
  #   figure <- figure + 
  #     coord_cartesian(default = TRUE, 
  #                     expand  = TRUE) +
  #     scale_y_continuous(#limits=c(-dy/5,NA),
  #       breaks = seq(0,1,by=0.25),
  #       limits = c(0, NA),
  #       expand = expansion(add = c(0.01, 0.1),
  #                          mult = mult),
  #       #limits = ifelse(percent, c(0, 100), NULL),
  #       labels = ifelse(percent, scales::percent, scales::number))
  # }
  
  
  return(figure)
  
}

plot_data <- function(input, 
                      x_summaries,
                      main_scenarios = NULL){
  
  dat <- x_summaries  %>%
    inner_join(input) %>% # should really carry these through when summarising
    mutate(second_test_delay_ = 
             ifelse(is.na(second_test_delay),
                    0,
                    second_test_delay),
           time_in_iso = 
             first_test_delay + 
             second_test_delay_+
             screening)
  
  
  if (!is.null(main_scenarios)){
    main_scenarios %<>% dplyr::select(-one_of("released_test")) %>% distinct
    dat <- left_join(dat, main_scenarios)
  }
  
  dat %>%  
    #tidyr::nest(data = -c(first_test_delay, second_test_delay)) %>%
    dplyr::mutate(delays = paste(first_test_delay, "&",
                                 first_test_delay + second_test_delay)) %>%
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
           #input,
           main_scenarios = NULL,
           plot = TRUE,
           log_scale = FALSE,
           #fixed = TRUE,
           text_size = 2.5,
           #trav_vol_manual = NULL,
           xlab = "Days in quarantine\n(including 1 day delay on testing results)",
           sum = FALSE,
           y_labels = NULL, # MUST BE PASSED IN!!!
           # pass in y_vars as a named list
           faceting = NULL,
           dir  = stringi::stri_rand_strings(1, 8),
           base = stringi::stri_rand_strings(1, 8)){
    
    if (!dir.exists(paste0("results/",dir))){
      dir.create(paste0("results/",dir))
    }
    
    #browser()
    all_grouping_vars <- all.vars(faceting)
    
    if (sum){
      y_labels <- sub("^Average", "Total", y_labels)
    } 
    
    #browser()
    
    ## need to do summaries over y_labels and x
    ## need a list length(y_labels) long, so perhaps an inner function mapping over?
    
    
    
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
          xlab      = "Days in quarantine\n(including 1 day delay on testing results)",
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