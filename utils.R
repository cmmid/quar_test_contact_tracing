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
  return(rlnorm(n, meanlog = meanlog, sdlog = sdlog))
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

calc_outcomes <- function(x, dat_gam){
  #browser()
  # generate required times for screening 
  # test 1: upon tracing
  #test 2: n days after 
  x <- mutate(x,
              first_test_t        = traced_t + first_test_delay,
              second_test_delay   = second_test_delay,
              second_test_t       = ifelse(traced_t > sec_exposed_t + second_test_delay,
                                           yes = traced_t,
                                           no = sec_exposed_t + second_test_delay)) %>% 
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
                                        sec_onset_t > traced_t &
                                        sec_onset_t < released_t ~
                                        "Symptomatic during quarantine",
                                      type == "symptomatic" & 
                                        sec_onset_t > sec_exposed_t &
                                        sec_onset_t < traced_t ~
                                        "Symptomatic before quarantine",
                                      TRUE ~ released_test),
           released_t = case_when(released_test == "Symptomatic during quarantine"~
                                    traced_t + pmax(sec_onset_t + post_symptom_window,
                                                    symp_end, max_mip),
                                  released_test =="Symptomatic before quarantine"~
                                    sec_exposed_t + pmax(sec_onset_t + post_symptom_window,
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
  

  incubation_times %<>% 
    mutate(
      onset     = exp_to_onset,
      symp_end  = ifelse(
        type == "asymptomatic",
        onset, # but really never matters because asymptomatics are never symptomatic!
        exp_to_onset + onset_to_recov),
      symp_dur  = symp_end - onset) # this doesn't get used because we know if someone's asymptomatic
  
  incubation_times %<>% gen_screening_draws
  
  return(incubation_times)
  
}


## just making sure the proportion of cases are secondary or not
make_sec_cases <- function(prop_asy, incubation_times){
  
  props <- c("asymptomatic" = prop_asy,
             "symptomatic"  = (1 - prop_asy))
  
  split_inc <- split(incubation_times,
                     incubation_times$type)
  
  res <- lapply(seq_along(props), function(x) sample_frac(split_inc[[x]], props[[x]]))
  
  do.call("rbind",res)
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
    
    #browser()
    
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
      rename("index_onset_t" = onset) %>% 
      mutate(index_testing_t    = index_onset_t + index_test_delay,
             index_result_t     = index_onset_t + index_test_delay + index_result_delay,
             traced_t           = index_onset_t + index_test_delay + index_result_delay +
               contact_info_delay + tracing_delay)
    
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
                     index_onset_t,
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
      ungroup() %>% 
      rename("sec_onset_t"=onset)
    
    ind_inc %<>%
      dplyr::select(-data)
    
    ind_inc %<>% 
      #rowwise %>%
      ## time of exposure of secondary cases is based on index's onset of symptoms
      ## it cannot be less than 0, hence the value of "a"
      ## it cannot be greater than some value... why?
      mutate(sec_exposed_t = index_onset_t - infect_shift + 
               rtgamma(n     = n(),
                       a     = pmax(0, infect_shift),
                       b     = infect_shift + index_testing_t, # for real?
                       shape = infect_shape,
                       rate  = infect_rate) 
      ) #%>% ungroup
    
    
    my_message("Shifting secondary cases' times relative to index cases' times")
    #exposure date relative to index cases exposure
    incubation_times_out <- ind_inc %>% 
      mutate(sec_onset_t     = sec_onset_t  + sec_exposed_t,
             symp_end_t  = symp_end  + sec_exposed_t)
    
    ## need to ditch dead columns
    rm(ind_inc)
    
    incubation_times_out <- left_join(input,
                                      incubation_times_out,
                                      by = c("index_test_delay", "delay_scaling"))
    
    #calc outcomes 
    my_message("Calculating outcomes for each traveller")
    incubation_times_out %<>% calc_outcomes(x       = .,
                                            dat_gam = dat_gam)
    
    my_message("Calculating when travellers released")
    incubation_times_out %<>% when_released(x = .)
    
    my_message("Transmission potential of released travellers")
    incubation_times_out %<>% transmission_potential
    
    return(incubation_times_out)
    
  }





transmission_potential <- function(x){
  
  x %<>% 
    mutate(
      q_exposed   = sec_exposed_t  - sec_onset_t + infect_shift,
      q_release   = released_t - sec_onset_t + infect_shift,
      q_traced    = traced_t   - sec_onset_t + infect_shift) 
  
  browser()
  
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
  
  # scale by infectivity_mass
  x %<>% 
    mutate_at(.vars = vars(infectivity_post,
                           infectivity_pre,
                           infectivity_quar),
              .funs = function(x,y){x/y}, y = .$infectivity_mass)
  
  x %<>% mutate(
    infectivity_total   =
      (infectivity_quar + infectivity_post + infectivity_pre),
    infectivity_averted = 1 - infectivity_total)
  
  return(x)
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

