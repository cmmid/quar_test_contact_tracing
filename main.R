
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")
source("lft_curves.R")

run_model <- function(
  n_sims          = 1000,
  n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
  n_ind_cases     = 10000,
  input,
  trajectories,
  seed            = 145,
  asymp_parms
){
  
  browser()
  
  set.seed(seed)
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))


sec_cases <- trajectories %>% 
  crossing(distinct(input)) %>%     
  mutate(index_testing_t  = index_onset_t + index_test_delay,
         sec_traced_t     = index_testing_t + test_to_tracing) %>% 
  mutate(sec_exposed_t = runif(n=n(),min = index_inf_start, max = pmin(index_inf_end, index_testing_t))) %>% 
  select(-data) %>% 
  mutate_at(.vars = vars(sec_onset_t),
            .funs = function(x,y){x + y}, y = .$sec_exposed_t) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso),
         adhering_symp=rbinom(n=n(),size = 1,prob = adherence_symp)) 

#browser()
# generate testing times
my_message("Calculating test times")
sec_cases %<>% 
  mutate(test_t = pmap(
    .l = list(
      multiple_tests = multiple_tests,
      tests          = tests,
      tracing_t      = sec_traced_t,
      sampling_freq  = sampling_freq,
      sec_exposed_t  = sec_exposed_t,
      quar_dur       = quar_dur,
      n_tests        = n_tests
    ),
    .f = test_times
  )) %>% 
  unnest(test_t) %>% 
  mutate(test_t = case_when(!is.na(assay) ~ test_t+lft_delivery_time,
                          TRUE            ~ test_t)) 


my_message("Generate possible missed days of testing")
sec_cases %<>% 
  unnest_longer(n_missed) %>% 
  group_by(ind_idx,sec_idx,n_missed,n_tests) %>% 
  nest() %>% 
  mutate(missed_days=pmap(.l=list(n_missed=n_missed,n_tests=n_tests),.f=missed_tests)) %>% 
  unnest(data) %>% 
  mutate(test_t = case_when(parse_number(test_no) %in% unlist(missed_days) ~ NA_real_,
                            TRUE                                           ~ test_t)) 

#calc outcomes 
my_message("Calculating outcomes for each secondary case")
sec_cases %<>% calc_outcomes(x  = .) %>% select(-c(m,rx,ry,u.x,u.y))

#find earliest positive test result
sec_cases %<>% 
  nest(test_t, test_q, test_p, test_no, test_label,  ct, screen) %>% 
  mutate(earliest_q      =  map(.f = earliest_pos, 
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
         symp_end_q = onset_q + post_symptom_window,
         test_iso_end_q = earliest_q + post_symptom_window
         )

# calculate remaining transmission potential averted by positive test
my_message("Calculating remaining transmission potential for each secondary case")
averted <- sec_cases %>% calc_overlap(.)
  
return(averted)

}

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `Daily testing` = 
      crossing(sampling_freq    = 1,
               tests            = TRUE,
               multiple_tests   = TRUE,
               n_tests          = c(3, 5, 7, 10), 
               lft_delivery_time  = c(0,2,4),
               assay            = c("Innova"
                                    #"Innova (+2.5 CT)",
                                    #"Innova (-2.5 CT)"
                                    ),
               quar_dur         = NA) %>% 
      mutate(n_missed          = list(c(NA,1,2))), 
    `Post-exposure quarantine only` = 
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA, 
               n_missed         = list(NA),
               quar_dur         = 10)#c(0, 5, 7, 10, 14))#,
    # `Post-exposure quarantine with LFA test` = 
    #   crossing(sampling_freq    = NA,
    #            tests            = TRUE,
    #            multiple_tests   = FALSE,
    #            n_tests          = NA,
    #            assay            = c("Innova"
    #                                 #"Innova (+2.5 CT)",
    #                                 #"Innova (-2.5 CT)"
    #                                ),
    #            quar_dur         = c(0, 5, 7, 10, 14)),
    # `Post-exposure quarantine with PCR test` = 
    #   crossing(sampling_freq    = NA,
    #            tests            = TRUE,
    #            multiple_tests   = FALSE,
    #            n_tests          = NA,
    #            assay            = "PCR", 
    #            quar_dur         = c(0, 5, 7, 10, 14))
  ) %>% 
    bind_rows(.id = "stringency")) %>% 
  mutate(post_symptom_window = 10,
           index_test_delay    = 1,  # time to entering quarantine (index cases)
           adherence_quar      = 1,# seq(0,1,by=0.25),
           adherence_iso       = 1,# seq(0,1,by=0.25),
           adherence_symp      = 1) %>% 
  crossing(test_to_tracing     = c(0,1,2)) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>%
  rowwise %>%
  group_split()

trajectories <- make_trajectories(n_ind = 1000,n_sec = 10,asymp_parms=asymp_fraction)

results_name <- "results_list"

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           trajectories= trajectories,
           seed = 1000
         )))

results_df <- get(results_name) %>% 
  bind_rows() %>%  
  select(-missed_days) %>% 
  as.data.frame() 

st=format(Sys.time(), "%Y%m%d")
write.fst(results_df,paste0("results/results_",st,"_main.fst"))


results_df %>% 
  filter(assay=="Innova",!is.na(n_tests),test_to_tracing==0,lft_delivery_time==0) %>% 
  mutate(n_missed=ifelse(is.na(n_missed),0,n_missed)) %>% 
  ungroup() %>% 
  mutate(test_number=extract_numeric(as.character(test_no))) %>% 
  mutate(test_number=replace_na(test_number,"None")) %>% 
  mutate(test_number=fct_relevel(test_number,"1","2","3","4","5","6","7","8","9","10","None"))%>%
  ggplot(aes(x=factor(n_tests),fill=fct_rev(as.factor(test_number))))+
  geom_bar(position="fill")+
  facet_grid(lft_delivery_time~n_missed+test_to_tracing,
             labeller = labeller(lft_delivery_time=function(x)paste(x,"days postage delay"),
                                 test_to_tracing=function(x)paste(x,"days tracing delay"),
                                 n_missed=function(x)paste(x,"days of tesing missed")))+
  scale_fill_manual(values = c("red",
                               viridis_pal(option = "viridis",begin=0.1,end=0.9,direction=-1)(10)))+
  labs(x="Days of testing",y="Proportion of infected individuals detected",fill="Detected by:")+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA))

ggsave("results/when_detected.png",height=120,width=210,dpi=400,units="mm")




