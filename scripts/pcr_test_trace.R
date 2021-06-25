
# Load required packages and utility scripts
source("scripts/packages.R")
source("scripts/utils.R")
source("scripts/plot_functions.R")
source("scripts/parameters.R")

run_model <- function(
  n_sims          = 1000,
  n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
  n_ind_cases     = 10000,
  input,
  trajectories,
  seed            = 145,
  asymp_parms
){
  

  #browser()
  
  set.seed(seed)
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))
  
  traj <- trajectories$traj %>% 
    select(-y) %>% 
    pivot_wider(names_from = name,values_from=x) %>% 
    select(-c(start,end))
  
  inf <- data.frame(prop_asy = rbeta(n = n_sims,
                                   shape1 = asymp_parms$shape1,
                                   shape2 = asymp_parms$shape2)) 



# Generate index cases' inc times
index_cases <- traj %>% 
  left_join(trajectories$models,by=c("idx","type")) %>% 
  filter(type=="symptomatic",!is.infinite(inf_start),!is.infinite(inf_end)) %>% 
  select(-c(data,m)) %>% 
  sample_n(n_sims) %>% 
  rename(ind_idx = idx) %>% 
  bind_cols(inf) %>% 
  #sample test result delay
  ## sample uniformly between 0 and 1 when 0.5...
  crossing(distinct(input, index_test_delay,test_to_tracing)) %>%     
  rename("index_onset_t"  = onset_t) %>% 
  mutate(index_testing_t  = index_onset_t + index_test_delay,
         sec_traced_t     = index_onset_t + index_test_delay + test_to_tracing) %>% 
  mutate(sec_exposed_t = runif(n=n(),min = inf_start,max=pmin(inf_end,index_testing_t)))




my_message("Generating secondary cases' trajectories and asymptomatic fraction")
sec_cases <- index_cases %>% 
  nest(data = -c(ind_idx,
                 prop_asy,
                 test_to_tracing,
                 index_onset_t,
                 index_test_delay,
                 index_testing_t,
                 sec_traced_t,
                 sec_exposed_t)) %>% 
  mutate(prop_asy    = as.list(prop_asy)) %>%
  mutate(sec_cases   = map(.x = prop_asy, 
                           .f  = ~make_sec_cases(as.numeric(.x),
                                                 traj,
                                                 n_sec_cases)
  ))

sec_cases %<>%
  unnest(prop_asy) %>%
  unnest(sec_cases) %>% 
  ungroup()%>% rename_at(.vars = vars(idx,onset_t),
                       .funs = ~paste0("sec_", .)) %>%
  dplyr::select(-data)


my_message("Shifting secondary cases' times relative to exposure to index")
#exposure date relative to index cases exposure
 sec_cases %<>% 
  mutate_at(.vars = vars(sec_onset_t),
            .funs = function(x,y){x + y}, y = .$sec_exposed_t) 


sec_cases <- left_join(input,
                       sec_cases,
                       by = c("index_test_delay", "test_to_tracing")) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_test=rbinom(n=n(),size = 1,prob = adherence_test),
         adhering_symp=rbinom(n=n(),size = 1,prob = adherence_symp)) 


# generate testing times
my_message("Calculating test times")
sec_cases %<>% 
  mutate(test_t = pmap(.l = list(multiple_tests=multiple_tests,
                                 tests=tests,
                                 tracing_t=sec_traced_t,
                                 sampling_freq=sampling_freq,
                                 sec_exposed_t=sec_exposed_t,
                                 quar_dur=quar_dur,
                                 max_tests=14,
                                 n_tests=n_tests),
                      .f = test_times)) %>% 
  unnest(test_t) 


#calc outcomes 
my_message("Calculating outcomes for each secondary case")
sec_cases %<>% calc_outcomes(x  = .) %>% select(-c(m,rx,ry,u.x,u.y))

#find earliest positive test result
sec_cases %<>% 
  nest(test_t, test_q, test_p, test_no, test_label, have_test, ct, detection_range, screen) %>% 
  mutate(earliest_q      =  map(.f = earliest_pos2, 
                                    .x = data)) %>% 
  unnest_wider(earliest_q) %>% 
  rename("earliest_q"=test_q)%>% 
  select(-data)

#shift other timings relative to onset
sec_cases %<>%
  mutate(exposed_q  = sec_exposed_t - sec_exposed_t,
         traced_q   = sec_traced_t - sec_exposed_t, 
         quar_end_q = pmax(traced_q,(exposed_q + quar_dur)),
         onset_q    = sec_onset_t - sec_exposed_t,
         symp_end_q = onset_q + iso_dur,
         test_iso_end_q = earliest_q + iso_dur,
         iso_release_test_q  = pmin(symp_end_q,test_iso_end_q)
         )

#release from isolation
sec_cases %<>% iso_release(x  = .) #%>% select(-c(m,rx,ry,u.x,u.y))

#browser()
# calculate remaining transmission potential averted by positive test
my_message("Calculating remaining transmission potential for each secondary case")
averted <- sec_cases %>% calc_overlap(.)
  
return(averted)

}

trajectories <- make_trajectories(n_cases = 10000)

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `1` = 
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA,
               test_exit_self_iso = FALSE,
               quar_dur         = 14,
               iso_dur          = 10),
    `2` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               test_exit_self_iso = FALSE,
               assay            = "LFA",
               sens_LFA         = "higher",
               quar_dur         = c(3,5,7,10),
               iso_dur          = 10),
    `3` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               test_exit_self_iso = TRUE,
               assay            = "LFA",
               sens_LFA         = "higher",
               quar_dur         = c(3,5,7,10),
               iso_dur          = c(3,5,7)),
    `4` = 
      crossing(sampling_freq    = 1,
               tests            = TRUE,
               multiple_tests   = TRUE,
               n_tests          = c(1,3,5,7,10), 
               test_exit_self_iso = FALSE,
               assay            = "LFA",
               sens_LFA         = "higher",
               quar_dur         = NA,
               iso_dur          = 10),
    `5` = 
      crossing(sampling_freq    = 1,
               tests            = TRUE,
               multiple_tests   = TRUE,
               n_tests          = c(1,3,5,7,10), 
               test_exit_self_iso = TRUE,
               assay            = "LFA",
               sens_LFA         = "higher",
               quar_dur         = NA,
               iso_dur          = c(3,5,7))
    
  ) %>% 
    bind_rows(.id = "stringency")) %>% 
  crossing(index_test_delay    = c(1),
           test_to_tracing       = c(0,3)) %>% 
  mutate(adherence_quar      = case_when(quar_dur==3~0.5,
                                         quar_dur==5~0.46,
                                         quar_dur==7~0.41,
                                         quar_dur==10~0.37,
                                         quar_dur==14~0.28),
         adherence_test  = case_when(iso_dur==3~1,
                                         iso_dur==5~0.97,
                                         iso_dur==7~0.93,
                                         iso_dur==10~0.86),
         adherence_symp = case_when(iso_dur==3~1,
                                         iso_dur==5~0.93,
                                         iso_dur==7~0.86,
                                         iso_dur==10~0.71)
         ) %>% 
  mutate(scenario=row_number())

input_split <-
  input %>%
  rowwise %>% 
  group_split()

results_name <- "results_list"


assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           trajectories=trajectories,
           seed = 1000,
           n_sec_cases = 10,
           n_sims = 500,
           asymp_parms = asymp_fraction
         )))

results_df <- get(results_name) %>% 
  bind_rows() %>% 
  as.data.frame() 

st=format(Sys.time(), "%Y%m%d_%H%M%S")
write.fst(results_df %>% select(-m),paste0("results/results_",st,"_act_a_1_and_10.fst"))
