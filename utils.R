my_message <- function(x, ...){
  message(paste(Sys.time(), x, sep = "    "), ...)
}


# is this not common to many scripts?
# move to utils.R
main_scenarios <-
  list(`low` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `moderate` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation")),
       `high` = 
         crossing(released_test = "Released after second test"),
       `maximum` = 
         crossing(released_test = c("Released after first test",
                                    "Released after mandatory isolation"))
  ) %>%
  bind_rows(.id = "stringency") %>%
  mutate(stage_released = "Infectious",
         stringency = fct_inorder(stringency)) 

probs        <- c(0.025,0.25,0.5,0.75,0.975)
lshtm_greens <- rev(c("#00BF6F","#0d5257"))

mv2gamma <- function(mean, var){
  list(shape = mean^2/var,
       rate  = mean/var,
       scale = var/mean) 
}

gamma2mv <- function(shape, rate=NULL, scale=NULL){
  if (is.null(rate)){
    rate <- 1/scale
  }
  
  list(mean = shape/rate,
       var  = shape/rate^2)
}


time_to_event <- function(n, mean, var){
  if (var > 0){
    parms <- mv2gamma(mean, var)
    return(rgamma(n, shape = parms$shape, rate = parms$rate))
  } else{
    return(rep(mean, n))
  }
}

time_to_event_lnorm <- function(n, meanlog, sdlog){
  return(rlnorm(n, meanlog = meanlog, sdlog = sdlog))
}




#Nishiura et al. 2020 serial interval
mean_si=4.7
sd_si=2.9

si <- distcrete::distcrete("lnorm",
                           meanlog=log(mean_si),
                           sdlog=log(sd_si),
                           interval = 1,
                           w = 0)

#He et al. Infectivity profile
infect_shape=97.18750 
infect_rate=3.71875
infect_shift=25.62500

# list of pathogens that may be worth considering as sensitivity


rriskDistributions::get.lnorm.par(q = c(5.1, 11.5),
                                  p = c(0.5, 0.975),
                                  plot = F,
                                  show.output = F) %>%
  as.list %>%
  set_names(., c("mu_inc", "sigma_inc")) -> inc_parms

pathogen <- list(
  symptomatic = 
    # review paper Byrne et al. (2020) https://doi.org/10.1101/2020.04.25.20079889
    # define T1 as infection to beginning of presymptomatic infectious period
    
    
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      
      # Li et al https://www.nejm.org/doi/full/10.1056/nejmoa2001316
      {c(9.1, 14.7)} %>% 
        set_names(., c("mu_inf", "sigma_inf"))),
  
  
  
  asymptomatic = 
    append(
      # https://www.acpjournals.org/doi/10.7326/M20-0504
      inc_parms,
      # https://doi.org/10.1101/2020.04.25.20079889
      list(
        mu_inf    =  6,
        sigma_inf = 12))) %>%
  map(~data.frame(.x), .id = "type")

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(q = c(0.24,  0.38),
                                                   p = c(0.025, 0.975), 
                                                   show.output = F, plot = F) %>%
  {list(shape1 = .[["shape1"]],
        shape2 = .[["shape2"]])}


gen_screening_draws <- function(x){
  #browser()
  
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x, dat_gam){
  #browser()
  # generate required times for screening 
  # test 1: upon tracing
  #test 2: n days after 
  x <- mutate(x,
              first_test_t        = traced_t + first_test_delay,
              second_test_delay   = second_test_delay,
              second_test_t       = ifelse(traced_t > exposed_t + second_test_delay,
                                           yes = traced_t,
                                           no = exposed_t + second_test_delay)) %>% 
    #if tests are on the same day, don't have the first test
    mutate(first_test_t = ifelse(second_test_t-first_test_t<1,
                                 yes = NA,
                                 no = first_test_t))
  
  # what's the probability of PCR detection at each test time?
  x <- mutate(x,first_test_p = 
                c(predict(object  = dat_gam,
                          type    = "response",
                          newdata = data.frame(day = first_test_t))),
              second_test_p = 
                c(predict(object  = dat_gam,
                          type    = "response",
                          newdata = data.frame(day = second_test_t)))
  ) 
  
  # asymptomatic infections have a lower detectability
  x <- mutate_at(x, 
                 .vars = vars(ends_with("test_p")),
                 .funs = ~ifelse(type == "asymptomatic",
                                 0.62 * .,
                                 .))
  # if pre-flight PCR is before infection, it's impossible to be detectable
  # this assumes 100% specificity
  # if we look at non-infected travellers in future this will change
  x <- mutate(x,
              first_test_p = ifelse(first_test_t < 0,
                                    0,
                                    first_test_p))
  
  # make comparisons of random draws to screening sensitivity
  # conduct symptomatic screening as well
  x <-
    mutate(x,
           first_test_label       = detector(pcr = first_test_p,  u = screen_1),
           second_test_label      = detector(pcr = second_test_p, u = screen_2))
}

when_released <- function(x){
  #browser()
  mutate(x, released_test = case_when(
    
    !screening    ~ 
      "Released after mandatory isolation",
    
    screening & is.na(first_test_label) & !second_test_label ~
      "Released after one test",
    
    screening & !first_test_label & !second_test_label     ~
      "Released after two tests",
    
    first_test_label                            ~
      "Released after one test + mandatory quarantine",
    
    !first_test_label  & second_test_label      ~
      "Released after two tests + mandatory quarantine",
    
    TRUE                                        ~ 
      "Mandatory quarantine"
  ),
  released_t = case_when(
    
    released_test == "Released after mandatory isolation"                   ~
      second_test_t, 
    
    released_test == "Released after one test"                           ~ 
      second_test_t + results_delay,
    
    released_test == "Released after two tests"                           ~ 
      second_test_t + results_delay,
    
    released_test == "Released after one test + mandatory quarantine"     ~ 
      second_test_t  + max_mip,
    
    released_test == "Released after two tests + mandatory quarantine"    ~
      second_test_t + max_mip,
    
    released_test == "Mandatory quarantine"                                 ~
      traced_t + max_mip)) %>% 
    mutate(released_test =  case_when(type == "symptomatic" & 
                                        onset > traced_t &
                                        onset < released_t ~
                                        "Symptomatic during quarantine",
                                      type == "symptomatic" & 
                                        onset > exposed_t &
                                        onset < traced_t ~
                                        "Symptomatic before quarantine",
                                      TRUE ~ released_test),
           released_t = case_when(released_test == "Symptomatic during quarantine"~
                                    traced_t + pmax(onset + post_symptom_window,
                                                    symp_end, max_mip),
                                  released_test =="Symptomatic before quarantine"~
                                    exposed_t + pmax(onset + post_symptom_window,
                                                     symp_end, max_mip),
                                  TRUE ~ released_t))
}

stage_when_released <- function(x){
  mutate(x, stage_released = as.factor(case_when(
    is.na(released_t)         ~ "Prevented from boarding",
    released_t < inf_end      ~ "Infectious",
    released_t >= inf_end     ~ "Post-infectious"
  ))) %>% 
    mutate(
      # how many days infectious after tracing
      days_released_inf = 
        case_when(
          stage_released == "Post-infectious" ~ 0,
          type           == "asymptomatic"    ~ inf_end - pmax(released_t, inf_start),
          type           == "symptomatic"     ~ onset   - pmax(released_t, inf_start)),
      # how many days infectious before tracing
      days_prior_inf =
        case_when(
          traced_t > inf_end   ~ inf_dur,
          traced_t < inf_start ~ 0,
          TRUE                 ~ traced_t - inf_start))
  
}


detector <- function(pcr, u = NULL, spec = 1){
  
  if (is.null(u)){
    u <- runif(n = length(pcr))
  }
  
  # true positive if the PCR exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- pcr > u
  
  # false positive if in the top (1-spec) proportion of random draws
  FP <- (pcr == 0)*(runif(n = length(pcr)) > spec)
  
  return(TP | FP)
}


make_plot_labels <- function(x){
  mutate(x, country = paste("Origin:",
                            ifelse(country == "United States of America",
                                   "USA", country))) %>%
    mutate_at(.vars = vars(contains("test_delay"),
                           pre_board_screening,
                           max_quarantine),
              .funs = function(x){ 
                sub(pattern = "1 days",
                    x = paste(x, "days"),
                    replacement = "1 day")
              })  %>%
    mutate(
      pre_board_screening = ifelse(pre_board_screening == "NA days",
                                   "No syndromic screening",
                                   paste("Syndromic screening", 
                                         pre_board_screening,
                                         "prior to departure"))) %>%
    mutate(max_quarantine = ifelse(max_quarantine == "0 days",
                                   "No isolation",
                                   paste("Isolate for",
                                         max_quarantine))) %>%
    mutate_at(.vars = vars(pre_board_screening,
                           matches("(first|second)_test_delay")),
              .funs = function(x){ifelse(grepl(pattern = "NA days",
                                               x = x),
                                         NA_character_,
                                         x)}) %>%
    mutate(stringency = stringr::str_to_title(stringency),
           stringency = paste(stringency, "stringency") ) %>%
    unite(col = "label",
          country,
          stringency,
          pre_board_screening,
          max_quarantine,
          sep = "\n") %>%
    dplyr::select(label) %>%
    {bind_cols(x, .)}
}


make_delay_label <- function(x,s){
  paste(na.omit(x), s)
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_incubation_times <- function(n_travellers,
                                  pathogen,
                                  asymp_parms){
  #browser()
  incubation_times <- crossing(i  = 1:n_travellers,
                               type = c("symptomatic",
                                        "asymptomatic") %>%
                                 factor(x = .,
                                        levels = .,
                                        ordered = T)) %>%
    mutate(idx=row_number()) %>% 
    dplyr::select(-i) %>% 
    split(.$type) %>%
    map2_df(.x = .,
            .y = pathogen,
            ~mutate(.x,
                    exp_to_onset   = time_to_event_lnorm(n = n(),
                                                         meanlog = .y$mu_inc, 
                                                         sdlog   = .y$sigma_inc),
                    onset_to_recov = time_to_event(n = n(),
                                                   mean = .y$mu_inf, 
                                                   var  = .y$sigma_inf))) 
  
  source("wolfel.R")
    
  
  incubation_times %<>% 
    mutate(
      onset     = exp_to_onset,
      symp_end  = ifelse(
        type == "asymptomatic",
        onset, # but really never matters because asymptomatics are never symptomatic!
        exp_to_onset + onset_to_recov),
      #inf_dur   = inf_end - inf_start,
      symp_dur  = symp_end - onset)
  
  incubation_times %<>% gen_screening_draws
  
  return(incubation_times)
  
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, incubation_times){
  #browser()
  
  #browser()
  props <- c("asymptomatic" = prop_asy,
             "symptomatic"  = (1 - prop_asy))
  
  split_inc <- split(incubation_times,
                     incubation_times$type)
  
  res <- lapply(seq_along(props), function(x) sample_frac(split_inc[[x]],props[[x]]))
  
  res <- do.call("rbind",res)
  res
}

make_arrival_scenarios <- function(input, 
                                   inf_arrivals, 
                                   incubation_times){
  #source('kucirka_fitting.R', local=T)
  
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  return(arrival_scenarios)
  
}

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


make_plots <- function(
  x, 
  input,
  main_scenarios = NULL,
  log_scale = FALSE,
  #fixed = TRUE,
  text_size = 2.5,
  #trav_vol_manual = NULL,
  xlab = "Days in quarantine\n(including 1 day delay on testing results)",
  sum = FALSE,
  y_var = "days_released_inf",
  faceting = NULL){
  
  #browser()
  
  ylabA = "Number of infectious persons\nreleased per index case"
  
  if (sum){
    ylabB = 
      "Total number of person-days infectiousness\nremaining for released secondary case"
  } else {
    ylabB = 
      "Number of days infectiousness\nremaining per released secondary case"
  }
  
  
  all_grouping_vars <- all.vars(faceting)
  
  # ... should be country, type, pre_board_screening, etc.
  x_summaries <- x %>%
    filter(stage_released == "Infectious") %>%
    make_released_quantiles(., all_grouping_vars) %>%
    inner_join(input)
  
  
  # why are the plo medians coming out backwars in Low?
  figA_data <- plot_data(input, 
                         x_summaries,
                         main_scenarios)
  
  
  # deal with the join and stage_released # need to have scenario guide the labelling
  
  # browser()
  figA <- figA_data %>% 
    #filter(pre_board_screening == "None") %>% 
    # # this filter should be done outside the function
    make_release_figure(
      x_summaries = .,
      input     = input,
      text_size = text_size,
      xlab      = xlab,
      ylab      = ylabA, 
      faceting  = faceting) 
  
  ## person-days
  
  x_days_summaries <- 
    make_released_time_quantiles(x, 
                                 y_var = y_var,
                                 all_grouping_vars,
                                 sum = sum)
  
  
  
  figB_data <- plot_data(input = input, 
                         x_summaries = 
                           x_days_summaries,
                         main_scenarios)
  
  figB <- figB_data %>% 
    #filter(pre_board_screening == "None") %>%  
    make_release_figure(
      x = .,
      input=input,
      xlab = "Days in quarantine\n(including 1 day delay on testing results)",
      text_size = text_size,
      ylab = ylabB,
      faceting = faceting) 
  
  
  fig <- figA + figB + plot_layout(ncol = 1, guide = "collect") +
    plot_annotation(tag_levels = "A",
                    theme = theme(legend.position = "bottom"))
  
  return(fig)
  
}



make_released_quantiles <- function(x, vars){
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  
  dots <- append(dots1, dots2)
  
  x_count <- x %>%
    dplyr::ungroup(.) %>%
    dplyr::select(!!! dots) %>%
    dplyr::group_by_all(.) %>%
    dplyr::count(.) 
  
  x_count %>%
    dplyr::ungroup(.) %>%
    dplyr::select(-n) %>%
    dplyr::ungroup %>%
    as.list %>%
    map(unique) %>%
    expand.grid %>%
    dplyr::left_join(x_count) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    tidyr::nest(data = c(sim, n)) %>%
    dplyr::mutate(
      Q = purrr::map(
        .x = data,
        .f = ~quantile(.x$n, probs = probs)),
      M = purrr::map_dbl(
        .x = data, 
        .f = ~mean(.x$n))) %>%
    tidyr::unnest_wider(Q) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup(.)
}

make_released_time_quantiles <- function(x, y_var, vars, sum = FALSE){
  #browser()
  
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  y_var <- as.name(y_var)
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, y_var) %>%
      group_by_at(.vars = vars(-y_var)) %>%
      summarise(y_var = sum(y_var, na.rm=T))
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, !! y_var) %>%
    dplyr::filter( !!y_var > 0)
  
  x_days %>%
    nest(data = c(!!y_var, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x[[y_var]],
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x[[y_var]]))) %>%
    unnest_wider(Q) %>%
    dplyr::select(-data)
  
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


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  return(gamma2mv(ans$estimate[["shape"]],
                  ans$estimate[["rate"]]))
}

run_analysis <- 
  function(n_sims          = 1000,
           n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
           n_ind_cases     = 10000,
           input,
           seed            = 145,
           P_c, P_r, P_t,
           dat_gam,
           asymp_parms){       # a list with shape parameters for a Beta
    
    browser()
    
    message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))
    
    #browser()
    set.seed(seed)
    
    my_message("Generating incubation times")
    
    #Generate incubation periods to sample
    incubation_times <- make_incubation_times(
      n_travellers = n_ind_cases,
      pathogen     = pathogen,
      asymp_parms  = asymp_parms)
    
    my_message("Generating asymptomatic fractions")
    inf <- data.frame(prop_asy = rbeta(n = n_sims,
                                       shape1 = asymp_parms$shape1,
                                       shape2 = asymp_parms$shape2)) 
    
    
    my_message("Generating index cases' transmissions")
    # Generate index cases' inc times
    ind_inc <- incubation_times %>% 
      filter(type=="symptomatic") %>% 
      sample_n(n_sims) %>% 
      mutate(sim = seq(1L, n_sims, by = 1L)) %>% 
      bind_cols(inf) %>% 
      #sample test result delay
      ## sample uniformly between 0 and 1 when 0.5...
      mutate(index_result_delay = time_to_event(n = n(),
                                                mean = P_r[["mean"]],
                                                var  = P_r[["var"]])) %>% 
      #sample contact info delay
      mutate(contact_info_delay = time_to_event(n = n(),
                                                mean = P_c[["mean"]],
                                                var  = P_c[["var"]])) %>% 
      #sample tracing delay
      mutate(tracing_delay      = time_to_event(n = n(),
                                                mean = P_t[["mean"]],
                                                var  = P_t[["var"]])) %>% 
      #add index test delay (assume 2 days post onset)
      #crossing(distinct(input, index_test_delay, delay_scaling, waning)) %>%     
      crossing(distinct(input, index_test_delay, delay_scaling, waning)) %>%     
      mutate_at(.vars = vars(tracing_delay, contact_info_delay, index_result_delay),
                .funs = ~(. * delay_scaling)) %>%
      rename("index_onset" = onset) %>% 
      mutate(index_testing_t    = index_onset + index_test_delay,
             index_result_t     = index_onset + index_test_delay + index_result_delay,
             traced_t           = index_onset + index_test_delay + index_result_delay +
               contact_info_delay + tracing_delay)
    
    rm(list = c("P_t", "P_r", "P_c", "inf"))
    
    my_message("Generating secondary cases' incubation times")
    
    #Generate secondary cases
    sec_cases <- make_incubation_times(
      n_travellers = n_sec_cases,
      pathogen     = pathogen,
      asymp_parms  = asymp_parms)
    
    #browser()
    my_message("Generating secondary cases' exposure times")
    ind_inc %<>% 
      nest(data = -c(sim,
                     prop_asy,
                     index_onset,
                     index_test_delay,
                     index_result_delay,
                     contact_info_delay,
                     tracing_delay,
                     index_testing_t,
                     traced_t,
                     delay_scaling))
    
    ind_inc %<>% 
      mutate(prop_asy    = as.list(prop_asy)) %>%
      mutate(sec_cases   = map(.x = prop_asy, 
                               .f  = ~make_sec_cases(as.numeric(.x),
                                                     sec_cases)
      ))
    
    rm(sec_cases)
    
    ind_inc %<>%
      unnest(prop_asy) %>%
      unnest(sec_cases) %>% 
      ungroup() 
    
    ind_inc %<>%
      dplyr::select(-data)
    
    ind_inc %<>% 
      #rowwise %>%
      ## time of exposure is based on index's onset of symptoms
      ## it cannot be less than 0, hence the value of a
      ## it cannot be greater than some value... why?
      mutate(exposed_t = index_onset - infect_shift + 
               rtgamma(n     = n(),
                       a     = pmax(0, infect_shift - index_onset),
                       b     = infect_shift + index_testing_t, # for real?
                       shape = infect_shape,
                       rate  = infect_rate) 
      ) #%>% ungroup
    
    
    my_message("Shifting secondary cases' times relative to index cases' times")
    #exposure date relative to index cases exposure
    incubation_times_out <- ind_inc %>% 
      mutate(onset     = onset     + exposed_t,
             symp_end  = symp_end  + exposed_t) 
    
    ## need to ditch dead columns
    rm(ind_inc)
    
    incubation_times_out <- left_join(input,
                                      incubation_times_out,
                                      by = c("index_test_delay", "delay_scaling"))
    
    
    #source('kucirka_fitting.R',local=T)  
    
    #calc outcomes 
    my_message("Calculating outcomes for each traveller")
    incubation_times_out %<>% calc_outcomes(x       = .,
                                            dat_gam = dat_gam)
    
    #when released
    my_message("Calculating when travellers released")
    incubation_times_out %<>% when_released(x = .)
    
    #stage of infection when released
    #my_message("Calculating infection status on release")
    # incubation_times_out %<>% stage_when_released()
    
    my_message("Transmission potential of released travellers")
    incubation_times_out %<>% transmission_potential
    
    return(incubation_times_out)
    
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

transmission_potential <- function(x){
  browser()
  
  x %<>% 
    mutate(
      onset_sec   = onset      - exposed_t, # index onset less exposure
      release_sec = released_t - exposed_t, 
      q_exposed   = exposed_t  - onset + infect_shift,
      q_release   = released_t - onset + infect_shift,
      q_traced    = traced_t   - onset + infect_shift) 
  
  x %<>%
    mutate(
      infectivity_mass        = pgamma(q     = q_exposed,
                                       shape = infect_shape,
                                       rate  = infect_rate,
                                       lower.tail = F),
      infectivity_post        = pgamma(q     = q_release, 
                                       shape = infect_shape,
                                       rate  = infect_rate, 
                                       lower.tail = F),
      infectivity_pre         = pgamma(q     = q_traced,
                                       shape = infect_shape,
                                       rate  = infect_rate),
    ) 
  
  #browser()
  
  if (all(x$waning == "waning_none")){
    # do pgamma
    x %<>% mutate(infectivity_quar = 0)
  } else {
    x %<>%
      mutate(
        infectivity_quar    =
          pmap_dbl(.l = list(q_traced,
                             q_release, 
                             waning),
                   .f = ~integrate(
                     f = function(x){
                       dgamma(x, shape = infect_shape, rate  = infect_rate) * 
                         (1 - get(..3)(x - ..1))
                     },
                     lower = ..1,
                     upper = ..2)$value))
  }
  
  # scale by infecivity_mass
  x %<>% 
    mutate_at(.vars = vars(infectivity_post,
                           infectivity_pre,
                           infectivity_quar),
              .funs = function(x,y){x/y}, y = .$infectivity_mass)
  
  x %<>% mutate(
    infectivity_total   =
      (infectivity_quar + infectivity_post + infectivity_pre),
    infectivity_averted = 1 - infectivity_total)
  
  # , 
  #               infectivity_pre = ptrunc(spec="gamma",
  #                                        a= traced_t -onset +infect_shift - 10,
  #                                        b= traced_t -onset +infect_shift + 10,
  #                                        q = traced_t - onset + infect_shift,
  #                                         shape = infect_shape,
  #                                         rate  = infect_rate,
  #                                         lower.tail = TRUE))
  
  return(x)
}

#function with if statements removed
ptrunc_v <- function (q, spec, a = -Inf, b = Inf, ...) 
{
  #browser()
  tt <- q
  aa <- rep(a, length(q))
  bb <- rep(b, length(q))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt <- G(apply(cbind(apply(cbind(q, bb), 1, min), aa), 1, 
                max), ...)
  tt <- tt - G(aa, ...)
  G.a <- G(aa, ...)
  G.b <- G(bb, ...)
  result <- tt/(G(bb, ...) - G(aa, ...))
  return(result)
}


rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  
  p_b <- pgamma(q = b, shape = shape, rate = rate)
  p_a <- pgamma(q = a, shape = shape, rate = rate)
  
  u   <- runif(n = n, min = p_a, max = p_b)
  q   <- qgamma(p = u, shape = shape, rate = rate)
  
  return(q)
}


check_unique_values <- function(df, vars){
  # given a data frame and a vector of variables to be used to facet or group, 
  # which ones have length < 1?
  
  l <- lapply(X = vars, 
              FUN =function(x){
                length(unique(df[, x]))
              })
  
  vars[l > 1]
  
}


waning_piecewise_linear <- function(x, ymax, ymin, k, xmax){
  
  if (ymin == ymax){
    Beta = c(0, ymin)
  } else {
    
    Beta <- solve(a = matrix(data = c(xmax, 1,
                                      k,    1),    ncol = 2, byrow = T),
                  b = matrix(data = c(ymin, ymax), ncol = 1))
  }
  
  (x >= 0)*pmin(ymax, pmax(0, Beta[2] + Beta[1]*x))
  
}

waning_points <- function(x, X, Y, log = FALSE){
  
  if (length(X) != length(Y)){
    stop("X and Y must be same length")
  }
  
  if (length(Y) == 1){
    return(rep(Y, length(x)))
  }
  
  if (log){
    Y <- log(Y)
  }
  
  Beta <- solve(a = cbind(X, 1), b = matrix(Y,ncol=1))
  
  Mu <- Beta[2] + Beta[1]*x
  if (log){
    Mu <- exp(Mu)
  }
  (x >= 0)*pmax(0, Mu)
  
}

