## functions used for ploting

# colours
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

released_palette <-
  c("Released after mandatory isolation"          = "#f2a19b",
    "Released after one negative test"            = "#2ac2db",
    "Released after two negative tests"           = "#5682c4",
    "Released after tests + mandatory quarantine" = "#df8df0")


#extrafont::loadfonts()
pdf.options(useDingbats=FALSE)


test_labeller <- function(x){
  mutate(x,
         stringency = factor(stringency,
                             levels = c("none",
                                        "one",
                                        "two"),
                             labels = c("No test",
                                        "One test",
                                        "Two tests"),
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


make_release_figure <- function(x,
                                input,
                                xlab = NULL,
                                text_size = 2.5,
                                text_angle = 45,
                                h_just = 0,
                                log_scale = FALSE,
                                y_var = "infectivity_averted",
                                hline = 0,
                                faceting = NULL,
                                percent = FALSE){
  
  x %<>% bind_rows %>% test_labeller  %>%
    mutate(released_test = factor(released_test,
                                  levels = names(released_palette),
                                  ordered = T)) # should this be in the facet call?
  
  x$y_var <- x[[y_var]]
  
  figure <- x %>%
    ggplot(data=., aes(x = y_var)) +
    geom_histogram(binwidth = 0.1, center = 0.05, 
                   aes(fill = released_test,
                       y = ..density..)
    ) +
    #facet_grid(
    facet_nested(nest_line = TRUE,
                 facets = faceting,
                 labeller = labeller(index_test_delay = index_test_labeller,
                                     delay_scaling    = delay_scaling_labeller,
                                     waning           = waning_labeller,
                                     stringency       = str_to_title)) +
    theme_minimal() +
    scale_fill_manual(values = released_palette,
                      name = "Release from quarantine") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +
    scale_x_continuous(limits = c(0,1), breaks = pretty_breaks(n=5))
  
  
  figure <- figure +
    theme(axis.ticks = element_line(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(fill=NA),
          legend.position = "bottom",
          strip.placement = "outside") +
    ggplot2::labs(x = xlab,
                  y = "Density") 
  
  
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
      
      faceting_lhs <- lhs(faceting)
      faceting_rhs <- rhs(faceting)
      
      faceting_new <- as.formula(
        paste(
          paste(formula.tools::get.vars(faceting_lhs), collapse = " + "),
          "~",
          paste(formula.tools::get.vars(faceting_rhs), collapse = " + "),
          "+",
          "second_test_delay",
          collapse = " ") %>% sub(pattern = "( \\. \\+ | \\+ \\. )",
                                  replacement = " ", x = .)
      )
      
      
      
      # avoid summarising the data in order to put it in
      
      figs <- map2(
        .x = y_labels,
        .y = names(y_labels),
        .f = ~make_release_figure(
          x         = x,
          text_size = text_size,
          y_var     = .y,
          xlab      = .x,
          faceting  = faceting_new,
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