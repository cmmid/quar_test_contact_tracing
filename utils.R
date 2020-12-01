my_message <- function(x, ...){
  message(paste(Sys.time(), x, sep = "    "), ...)
}

# is this not common to many scripts?
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
  rlnorm(n, meanlog = meanlog, sdlog = sdlog)
}

gen_screening_draws <- function(x){
  n <- nrow(x)
  
  # generate screening random draws for comparison
  x <- mutate(x, 
              screen_1 = runif(n, 0, 1),  # on arrival
              screen_2 = runif(n, 0, 1))  # follow-up
}

# given infection histories above, what proportion of travellers end up being 
# caught at each step in the screening process?

calc_outcomes <- function(x){
  #browser()
  # generate required times for screening 

  x <- x %>% 
    mutate(test_t=ifelse(test_t<sec_traced_t,
                         sec_traced_t,
                         test_t))
  
  
  # what's the probability of detection at each test time given a value of CT?
  x_ <- x %>%
   inner_join(trajectories$models, by=c("sec_idx"="idx","type")) %>% 
    select(-data) %>% 
    mutate(test_q=test_t-sec_exposed_t,
           ct = map2_dbl(.f = predict,
                         .x = m,
                         .y = test_q)) %>% 
    mutate(detection_range = cut(ct,breaks = c(-Inf,27,30,35,Inf))) %>% 
    mutate(test_p=case_when(assay=="PCR"&detection_range=="(-Inf,27]" ~ 1,
                                    assay=="PCR"&detection_range=="(27,30]"   ~ 0.75,
                                    assay=="PCR"&detection_range=="(30,35]"   ~ 0.5,
                                    assay=="PCR"&detection_range=="(35, Inf]" ~ 0,
                                    assay=="LFA"&detection_range=="(-Inf,27]" ~ 0.9,
                                    assay=="LFA"&detection_range=="(27,30]"   ~ 0.5,
                                    assay=="LFA"&detection_range=="(30,35]"   ~ 0,
                                    assay=="LFA"&detection_range=="(35, Inf]" ~ 0)) %>% 
    mutate(screen      = runif(n(), 0, 1)) %>% 
    mutate(test_label  = detector(pcr = test_p,  u = screen))
    
    
    # # can't return a test prior to exposure
    # x_ <- mutate(x_,
    #              test_p = ifelse(test_t < sec_exposed_t,
    #                              NA, # this may need to be NA
    #                              test_p))
    # 
 
  return(x_)
}

when_released <- function(x){
  # browser()
  # NOT REVIEWED YET
  mutate(x, 
         released_test = case_when(
           
           stringency == "none" ~
             "Released after mandatory quarantine",
           
           is.na(first_test_label) & is.na(second_test_label) ~
             "Released after mandatory quarantine",
           
           is.na(first_test_label) & !second_test_label ~
             "Released after negative end of quarantine test",
           
           !first_test_label & !second_test_label       ~
             "Released after two negative tests",
           
           
           first_test_label | !second_test_label ~
             "Released after positive first test + mandatory isolation",
           
           !first_test_label & second_test_label | is.na(first_test_label) & second_test_label ~
             "Released after positive second test + mandatory isolation",
           
           TRUE                                         ~ 
             "ILLEGAL CONFIGURATION. Cannot have false first test and NA second test"
         ),
         released_t = case_when(
           
           released_test == "Released after mandatory quarantine"     ~
             second_test_t, 
           
           released_test == "Released after negative end of quarantine test"       ~ 
             second_test_t + results_delay * delay_scaling,
           
           released_test == "Released after two negative tests"     ~ 
             second_test_t + results_delay * delay_scaling,
           
           released_test == "Released after positive first test + mandatory isolation"     ~ 
             first_test_t + post_symptom_window,
           
           released_test == "Released after positive second test + mandatory isolation"     ~ 
             second_test_t + post_symptom_window)) %>% 
    mutate(released_test_symptomatic = 
             case_when(type == "symptomatic" & 
                         sec_onset_t >= index_traced_t &                         sec_onset_t < released_t ~
                         "Symptomatic during quarantine",
                       type == "symptomatic" & 
                         sec_onset_t < index_traced_t ~
                         "Symptomatic before quarantine",
                       type == "symptomatic" &
                         sec_onset_t >= released_t ~
                         "Symptomatic after quarantine",
                       TRUE ~ "Never symptomatic"))#,
  # released_t    = case_when(
  #   released_test_symptomatic == "Symptomatic during quarantine"~
  #     pmax(sec_onset_t + post_symptom_window,
  #          sec_symp_end_t),
  #   released_test_symptomatic == "Symptomatic before quarantine"~
  #     pmax(sec_onset_t + post_symptom_window,
  #          sec_symp_end_t,
  #          second_test_t),
  #   TRUE ~ released_t))
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
  
  TP | FP
}


make_delay_label <- function(x,s){
  paste(na.omit(x), s)
}


capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

make_trajectories <- function(n_cases){
  #simulate CT trajectories
  #browser()
  traj <- data.frame(idx=1:n_cases) %>% 
    crossing(start=0,type=c("symptomatic","asymptomatic") %>% 
               factor(x = .,
                      levels = .,
                      ordered = T)) %>% 
    mutate(end=case_when(type=="symptomatic" ~ rnorm(n=n(),mean=17,sd=2),
                         type=="asymptomatic" ~ rnorm(n=n(),mean=17*0.6,sd=1))) %>% 
    mutate(onset_t=rlnormTrunc(n=n(),
                            meanlog=1.63,
                            sdlog=0.5,
                            min=0,
                            max=end)) %>% 
    pivot_longer(cols=-c(idx,type),values_to = "x") %>% 
    mutate(y=case_when(name=="start"~40,
                       name=="end"~40,
                       name=="onset_t"~rnorm(n=n(),mean=25,sd=5))) 
  
  models <- traj %>%
    nest(-c(idx,type)) %>%  
    dplyr::mutate(
      # Perform loess calculation on each individual 
      m = purrr::map(data, loess,
                     formula = y ~ x),
      inf_period=purrr::map(.f=infectious_period,
                            .x=m)) %>% 
    unnest_wider(inf_period)
  
  return(list(traj=traj,models=models))
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, trajectories,n_sec_cases){

  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  res <- lapply(names(props), 
                function(x){
                  filter(trajectories, type == x) %>%
                    sample_frac(., size = props[[x]])
                })
  
  res <- do.call("rbind",res) %>% 
    sample_n(n_sec_cases)
  
}

make_arrival_scenarios <- function(input, 
                                   inf_arrivals, 
                                   incubation_times){
  #source('kucirka_fitting.R', local=T)
  
  arrival_scenarios <- crossing(input, inf_arrivals)
  
  # calculate outcomes of screening
  arrival_scenarios %<>% calc_outcomes(., dat_gam)
  
  arrival_scenarios
  
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

make_quantiles <- function(x, y_var,vars, sum = TRUE){
  #browser()
  dots1 <- rlang::exprs(sim, scenario)
  dots2 <- lapply(vars, as.name)
  y_var <- as.name(y_var)
  dots  <- append(dots1, dots2)
  
  if (sum){
    x <- x %>%
      dplyr::select(!!! dots, y_var) %>%
      group_by_at(.vars = vars(-y_var)) %>% 
      summarise(!!y_var := sum(!!y_var, na.rm=T),n=n()) %>% 
      mutate(!!y_var:=!!y_var/n)
  }
  
  x_days <- x %>%
    dplyr::select(!!! dots, !! y_var) 
  
  x_days %>%
    nest(data = c(!!y_var, sim)) %>%
    mutate(Q = purrr::map(.x = data, ~quantile( .x[[y_var]],
                                                probs = probs)),
           M = map_dbl(.x = data, ~mean(.x[[y_var]]))) %>%
    unnest_wider(Q) %>% 
    # mutate(Q = purrr::map(.x = data, ~bayestestR::hdi( .x[[y_var]],
    #                                             ci = c(0.95,0.5))),
    #        Mean = map_dbl(.x = data, ~mean(.x[[y_var]])),
    #        MAP = map_dbl(.x=data, ~as.numeric(bayestestR::point_estimate(.x[[y_var]],centrality="MAP")))) %>% 
    #unnest(Q) %>% 
    #pivot_wider(names_from = CI,values_from=c(CI_low,CI_high)) %>% 
    dplyr::select(-data) 
  
}


delay_to_gamma <- function(x){
  ans <- dplyr::transmute(x, left = t - 0.5, right = t + 0.5) %>%
    dplyr::mutate(right = ifelse(right == max(right), Inf, right)) %>%
    {fitdistcens(censdata = data.frame(.),
                 distr = "gamma", 
                 start = list(shape = 1, rate = 1))} 
  
  gamma2mv(ans$estimate[["shape"]],
           ans$estimate[["rate"]])
}

run_analysis <- 
  function(n_sims          = 1000,
           n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
           n_ind_cases     = 10000,
           input,
           seed            = 145,
           P_c, P_r, P_t,
           dat_gam,
           asymp_parms,
           return_full = TRUE,
           faceting = stringency ~ type){       # a list with shape parameters for a Beta
    
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
      #sample tracing delay (per sec case)
      mutate(tracing_delay      = time_to_event(n = n(),
                                                mean = P_t[["mean"]],
                                                var  = P_t[["var"]])) %>% 
      #add index test delay (assume 2 days post onset)
      crossing(distinct(input, index_test_delay, delay_scaling,test_to_tracing)) %>%     
      mutate_at(.vars = vars(index_result_delay,contact_info_delay,tracing_delay),
                .funs = ~(. * delay_scaling)) %>%
      rename("index_onset_t" = onset) %>% 
      mutate(index_testing_t    = index_onset_t + index_test_delay,
             index_result_t     = index_onset_t + index_test_delay + index_result_delay,
             index_traced_t     = index_onset_t + index_test_delay + test_to_tracing)#+ index_test_delay + index_result_delay +
               #contact_info_delay + tracing_delay)
    
    
    #rm(list = c("P_t", "P_r", "P_c", "inf"))
    
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
                     test_to_tracing,
                     index_onset_t,
                     index_test_delay,
                     index_result_delay,
                     contact_info_delay,
                     tracing_delay,
                     index_testing_t,
                     index_traced_t,
                     delay_scaling))
    
    ind_inc %<>% 
      mutate(prop_asy    = as.list(prop_asy)) %>%
      mutate(sec_cases   = map(.x = prop_asy, 
                               .f  = ~make_sec_cases(as.numeric(.x),
                                                     sec_cases)
      ))
    
    #rm(sec_cases)
    
    
    ind_inc %<>%
      unnest(prop_asy) %>%
      unnest(sec_cases) %>% 
      ungroup()
    
    ind_inc %<>% rename_at(.vars = vars(onset, symp_end, symp_dur,
                                        exp_to_onset, onset_to_recov),
                           .funs = ~paste0("sec_", .))
    
    ind_inc %<>%
      dplyr::select(-data)
    
    ind_inc %<>% 
      #rowwise %>%
      ## time of exposure of secondary cases is based on index's onset of symptoms
      ## it cannot be less than 0, hence the value of "a"
      ## it cannot be greater than some value... why?
      mutate(sec_exposed_t = index_onset_t - infect_shift + 
               rtgamma(n     = n(),
                       a     = infect_shift,
                       b     = infect_shift + index_testing_t - index_onset_t, 
                       shape = infect_shape,
                       rate  = infect_rate) 
      ) #%>% ungroup
    
    my_message("Shifting secondary cases' times relative to index cases' times")
    #exposure date relative to index cases exposure
    incubation_times_out <- ind_inc %>% 
      rename_at(.vars = vars(sec_onset, sec_symp_end),
                .funs = ~paste0(., "_t")) %>%
      mutate_at(.vars = vars(sec_onset_t, sec_symp_end_t),
                .funs = function(x,y){x + y}, y = .$sec_exposed_t) 
    

    
    ## need to ditch dead columns
    #rm(ind_inc)
    
     incubation_times_out <- left_join(input,
                                       incubation_times_out,
                                       by = c("index_test_delay", "delay_scaling","test_to_tracing")) %>% 
       mutate(adhering=rbinom(n=n(),size = 1,prob = adherence)) 
     
     incubation_times_out %<>%  
       mutate(test_t = map(.x = sampling_freq,
                           .f = test_times)) %>% 
       unnest(test_t)
    
    #calc outcomes 
    my_message("Calculating outcomes for each secondary case")
    incubation_times_out %<>% calc_outcomes(x       = .,
                                            dat_gam = dat_gam)
    
    my_message("Calculating when secondary cases released")
    incubation_times_out %<>% when_released
    
    my_message("Transmission potential of released secondary case")
    incubation_times_out %<>% transmission_potential
    
    #browser()
    
    if (return_full){
      my_message("Returning simulation results")
      return(incubation_times_out)
    } else {
      my_message("Calculating and returning simulation summary statistics")
      # pull into a function
      
      return(summarise_simulation(incubation_times_out,faceting))
      
    }
    
    
    
  }

transmission_potential <- function(x){
  #browser()
  x %<>% 
    mutate(
      q_exposed   = sec_exposed_t  - sec_onset_t + infect_shift,
      q_release   = released_t     - sec_onset_t + infect_shift,
      q_traced    = index_traced_t - sec_onset_t + infect_shift,
      q_onset     = sec_onset_t    - sec_onset_t + infect_shift,
      q_symp_end  = sec_symp_end_t - sec_onset_t + infect_shift,
      q_leave_iso = pmax(q_symp_end, q_onset + post_symptom_window,q_release))  
  
  
  #if you're symptomatic before quarantine, you enter self-isolation for 10 days
  
  x %<>%
    # mutate(trans_pot_post_exp = 
    #          pgamma(q     = q_exposed, 
    #                 shape = infect_shape,
    #                 rate  = infect_rate,
    #                 lower.tail = FALSE)) %>%
    # mutate(trans_pot_start_quar= pgamma(q     = q_traced, 
    #                                    shape = infect_shape,
    #                                    rate  = infect_rate,
    #                                    lower.tail = FALSE)) %>% 
    # mutate(trans_pot_end_quar =  pgamma(q     = q_release, 
    #                                     shape = infect_shape,
    #                                     rate  = infect_rate,
    #                                     lower.tail = FALSE)) %>% 
    # mutate(trans_pot_in_quar =   trans_pot_start_quar-trans_pot_end_quar) %>% 
    # mutate(trans_pot_start_symp = case_when(
    #   type=="symptomatic"~pgamma(q = onset_q, 
    #                              shape = infect_shape,
    #                              rate  = infect_rate,
    #                              lower.tail = FALSE),
    #   type=="asymptomatic"~NA_real_)) %>% 
    # mutate(trans_pot_end_symp   =  pgamma(q     = pmax(q_symp_end, q_onset + post_symptom_window), 
    #                                       shape = infect_shape,
    #                                       rate  = infect_rate,
    #                                       lower.tail = FALSE)) %>% 
    # mutate(trans_pot_in_symp = trans_pot_start_symp-trans_pot_end_symp) %>% 
    # mutate(trans_pot_pos_test = pgamma(q = earliest_q, 
    #                                    shape = infect_shape,
    #                                    rate  = infect_rate,
    #                                    lower.tail = FALSE)) %>% 
    # mutate_at(.vars = vars(trans_pot_pos_test,
    #                        trans_pot_post_symp),
    #           .funs = function(x,y){x/y},
    #           y = .$trans_pot_post_exp) %>% 
    # mutate(trans_pot_post_symp=
    #          case_when(
    #            type == "symptomatic" & adhering==1 ~ 
    #              trans_pot_post_symp,
    #            type == "symptomatic" & adhering==0 ~ 
    #              0
    #          )) %>% 
    # mutate(trans_pot_averted = case_when(
    #   type == "symptomatic" & (trans_pot_pos_test > trans_pot_post_symp) ~ 
    #     trans_pot_pos_test - trans_pot_post_symp,
    #   type == "symptomatic"  ~ 0,
    #   TRUE ~ trans_pot_pos_test)) # asymptomatics
  mutate(infectivity_pre=case_when(released_test_symptomatic == "Symptomatic before quarantine" ~
                                       pmap_dbl(.l = list(q_exposed,
                                                          q_onset),
                                                .f = ~integrate(
                                                  f = function(x){
                                                    dgamma(x, shape = infect_shape, rate = infect_rate)
                                                  },
                                                  lower = ..1,
                                                  upper = ..2)$value),
                                     released_test_symptomatic != "Symptomatic before quarantine" ~
                                       pmap_dbl(.l = list(q_exposed,
                                                          q_traced),
                                                .f = ~integrate(
                                                  f = function(x){
                                                    dgamma(x, shape = infect_shape, rate  = infect_rate)
                                                  },
                                                  lower = ..1,
                                                  upper = ..2)$value)),
           infectivity_quar=case_when(released_test_symptomatic == "Symptomatic before quarantine" ~
                                        pmap_dbl(.l = list(q_leave_iso,
                                                           q_release,
                                                           waning),
                                                 .f = ~integrate(
                                                   f = function(x){
                                                     dgamma(x, shape = infect_shape, rate  = infect_rate) *(get(..3)(x-..1))
                                                   },
                                                   lower = ..1,
                                                   upper = ..2)$value),
                                      released_test_symptomatic == "Symptomatic during quarantine" ~
                                        pmap_dbl(.l = list(q_traced,
                                                           q_onset,
                                                           waning),
                                                 .f = ~integrate(
                                                   f = function(x){
                                                     dgamma(x, shape = infect_shape, rate  = infect_rate) *(get(..3)(x-..1))
                                                   },
                                                   lower = ..1,
                                                   upper = ..2)$value),
                                      released_test_symptomatic %in% c("Symptomatic after quarantine") ~
                                        pmap_dbl(.l = list(q_traced,
                                                           q_release,
                                                           waning),
                                                 .f = ~integrate(
                                                   f = function(x){
                                                     dgamma(x, shape = infect_shape, rate  = infect_rate) *(get(..3)(x-..1))
                                                   },
                                                   lower = ..1,
                                                   upper = ..2)$value),
                                      released_test_symptomatic %in% c("Never symptomatic") ~
                                        pmap_dbl(.l = list(q_traced,
                                                           q_release,
                                                           adhering),
                                                 .f = ~integrate(
                                                   f = function(x){
                                                     dgamma(x, shape = infect_shape, rate  = infect_rate) 
                                                   },
                                                   lower = ..1,
                                                   upper = ..2)$value)),
           infectivity_iso=case_when(released_test_symptomatic != "Never symptomatic" ~ 
                                       pmap_dbl(.l = list(q_onset,
                                                          q_leave_iso,
                                                          adhering),
                                                .f = ~integrate(
                                                  f = function(x){
                                                    dgamma(x, shape = infect_shape, rate  = infect_rate)
                                                  },
                                                  lower = ..1,
                                                  upper = ..2)$value*..3),
                                     released_test_symptomatic == "Never symptomatic" ~ 
                                       0),
           infectivity_post = case_when(released_test_symptomatic == "Never symptomatic" ~ 
                                          pmap_dbl(.l = list(q_release),
                                                   .f = ~integrate(
                                                     f = function(x){
                                                       dgamma(x, shape = infect_shape, rate  = infect_rate)
                                                     },
                                                     lower = ..1,
                                                     upper = Inf)$value),
                                        released_test_symptomatic == "Symptomatic after quarantine" ~ 
                                          pmap_dbl(.l = list(q_release,
                                                             q_onset),
                                                   .f = ~integrate(
                                                     f = function(x){
                                                       dgamma(x, shape = infect_shape, rate  = infect_rate)
                                                     },
                                                     lower = ..1,
                                                     upper = ..2)$value),
                                        TRUE ~ 
                                          pmap_dbl(.l = list(q_leave_iso),
                                                   .f = ~integrate(
                                                     f = function(x){
                                                       dgamma(x, shape = infect_shape, rate  = infect_rate)
                                                     },
                                                     lower = ..1,
                                                     upper = Inf)$value)
           ),
         infectivity_post_exp =pgamma(q=q_exposed,
                                      shape=infect_shape,
                                      rate=infect_rate,
                                      lower.tail = F)
    ) %>% 
    mutate(infectivity_quar = ifelse(infectivity_quar<0, 0, infectivity_quar)) %>% 
    # mutate_at(.vars = vars(infectivity_pre,
    #                        infectivity_post,
    #                        infectivity_quar,
    #                        infectivity_iso),
    #           .funs = function(x,y){x/y},
    #           y = .$infectivity_post_exp) %>% 
    mutate(
      infectivity_spent   = infectivity_pre  + infectivity_post,
      infectivity_averted = infectivity_quar + infectivity_iso)
  
  return(x)
}



rtgamma <- function(n = 1, a = 0, b = Inf, shape, rate = 1, scale = 1/rate){
  #browser()
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



summarise_simulation <- function(x, faceting, y_labels = NULL){
  
  if(is.null(y_labels)){
    # if none specified, use all.
    y_labels_names <- grep(x=names(x), pattern="^trans_pot_", value = T)
  } else {
    y_labels_names <- names(y_labels)
  }
  
  all_grouping_vars <- all.vars(faceting)
  
  # if (!any(grepl(pattern = "type", x = all_grouping_vars))){
  #   all_grouping_vars <- c(all_grouping_vars, "type")
  # }
  
  x_summaries <-
    as.list(y_labels_names) %>%
    set_names(., .) %>%
    lapply(X = ., 
           FUN = function(y){
             make_quantiles(x,
                            y_var = y, 
                            vars = all_grouping_vars)})
  
  if (any(grepl(pattern = "type", x = all_grouping_vars))){
    
    x_summaries_all <- as.list(y_labels_names) %>%
      set_names(., .) %>%
      lapply(X = ., 
             FUN = function(y){
               make_quantiles(
                 mutate(x,
                        type = "all"),
                 y_var = y, 
                 vars = all_grouping_vars)})
    
    x_summaries <- map2(.x = x_summaries,
                        .y = x_summaries_all,
                        .f = ~bind_rows(.x, .y))
    
  }
  
  bind_rows(x_summaries, .id = "yvar")
  
}

calc_sensitivity <- function(model, x){
  browser()
  s <- predict(model,x)
  
  return(s)
}



read_results <- function(results_path){
  #browser()
  list(here::here("results", results_path, "results.rds"),
       here::here("results", results_path, "input.rds")) %>%
    map(read_rds) %>%
    map(bind_rows) %>%
    {inner_join(.[[1]], .[[2]])}
}


test_times <- function(tracing_t,sampling_freq = 7, max_time = 30, max_tests = NA){
 # browser()
  
  test1 <- tracing_t
  
  if (!is.na(max_tests)){
    max_time <- sampling_freq * (max_tests-1) + test1
  }
  
   if (is.na(sampling_freq)){
    test <- max_time
   } else{
     test  <- seq(from = test1, to = max_time, by=sampling_freq)
   }
  
  tests <- as.data.frame(test) %>% 
    rename("test_t" = test) %>% 
    mutate(test_no = paste0("test_", row_number()))
  
  return(tests)
}

earliest_pos <- function(df){
  #browser()
  x_q <- unlist(df[(df$test_label),"test_q"])
  
  if (length(x_q) == 0L){
    return(Inf)
  } else {
    return(min(x_q))
  }
}

earliest_pos2 <- function(df){
  #browser()
  x_q <- df %>% filter(test_label==TRUE) 
  
  if (nrow(x_q) == 0L){
    return(data.frame(test_no="None",test_q=Inf))
  } else {
    return((x_q %>% select(test_no,test_q) %>% slice_min(test_q)))
  }
}


infectious_period <- function(model){
  newdata <- data.frame(x=seq(from=0,to=30,0.1))
  
  newdata$y_pred <- predict(model,newdata)
  
  newdata %>% 
    mutate(infectious=y_pred<=30) %>% 
    filter(infectious) %>% 
    summarise(inf_start=min(x),
           inf_end=max(x))
}

calc_overlap <- function(a_min,a_max,b_min,b_max){
  Overlap(x=c(a_min,a_max),y=c(b_min,b_max))
}
