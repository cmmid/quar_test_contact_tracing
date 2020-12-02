
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("kucirka_fitting.R")
source("parameters.R")

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
  crossing(distinct(input, index_test_delay, delay_scaling,test_to_tracing)) %>%     
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
                 sec_exposed_t,
                 delay_scaling)) %>% 
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
                       by = c("index_test_delay", "delay_scaling","test_to_tracing")) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso)) 

#browser()
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
sec_cases %<>% calc_outcomes(x  = .)

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
         symp_end_q = onset_q + post_symptom_window,
         test_iso_end_q = earliest_q + post_symptom_window
         )
#browser()
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
               n_tests          = default_testing, 
               assay            = "LFA",
               sens_scaling     = c(0.67,1),
               quar_dur         = NA), 
    `Post-exposure quarantine only` = 
      crossing(sampling_freq    = NA,
               tests            = FALSE,
               multiple_tests   = FALSE,
               assay            = NA,
               n_tests          = NA,
               quar_dur         = c(0, default_testing[-1])),
    `Post-exposure quarantine with LFA test` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               assay            = "LFA",
               sens_scaling     = c(0.67,1),
               quar_dur         = c(0, default_testing[-1])),
    `Post-exposure quarantine with PCR test` = 
      crossing(sampling_freq    = NA,
               tests            = TRUE,
               multiple_tests   = FALSE,
               n_tests          = NA,
               sens_scaling     = 1,
               assay            = "PCR",
               quar_dur         = c(0, default_testing[-1]))
  ) %>% 
    bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window = 10,
           index_test_delay    = c(1),  # time to entering quarantine (index cases)
           delay_scaling       = c(1, 0.5, 0),
           adherence_quar      = c(0, 0.5,  1),
           adherence_iso       = c(0, 0.67, 1)) %>% 
  mutate(test_to_tracing       = 3*delay_scaling) %>% 
  filter(#adherence_iso==0.67,adherence_quar==0.5#,assay=="PCR"
         ) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>%
  rowwise %>%
  group_split()

trajectories <- make_trajectories(n_cases = 1000)

results_name <- "results_list"

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           trajectories=trajectories,
           seed = 1000,
           n_sec_cases = 10,
           n_sims = 100,
           asymp_parms = asymp_fraction
         )))

saveRDS(get(results_name),"results_20201202_all.rds")

assign(x=results_name,value=read_rds("results_20201201_test.rds"))

col_pal <- RColorBrewer::brewer.pal(n=4,name = "Dark2")

plot_a <- get(results_name)%>% 
  bind_rows() %>% 
  filter(#adherence_iso==1,
         #adherence_quar==1,
         delay_scaling==1,
         # assay!="PCR",
         !multiple_tests
  ) %>%
mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                          tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                          tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                          !tests~"Post-exposure quarantine only"
)) %>%
  group_by(ind_idx,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
            ) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_b <-get(results_name) %>% 
  bind_rows() %>% 
  filter(#adherence_iso==1,
         #adherence_quar==1,
         delay_scaling==1,
         #assay!="PCR",
         multiple_tests) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
            ) %>% 
  group_by(strategy,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  plotting_theme

plot_a+plot_b+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

  save_plot(dpi = 400, 
          device = "png",
          prefix = "daily_vs_end_quar_n_tests",
          base = "plot", 
          width = 300, 
          height = 150)


#In text
get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==1,
    adherence_quar==1,
    delay_scaling==1
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(sim,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(n=n(),
            prop=sum(trans_pot_averted)/n) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(iqr=paste0("(",iqr,")"),
         ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,strategy,quar_dur,`50%`,iqr,ui) %>% 
  htmlTable()

get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(sim,strategy,adherence_quar,assay,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(iqr=paste0("(",iqr,")"),
         ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,strategy,quar_dur,`50%`,iqr,ui) %>% 
  htmlTable()


### Sensitivity analyses ----
plot_a_delays <- get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    #delay_scaling==1,
    !multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_b_delays <- get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    #delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,assay,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,assay,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme

plot_a_delays+plot_b_delays+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "daily_vs_end_quar_delays",
          base = "plot", 
          width = 350, 
          height = 150)

get(results_name) %>% 
  bind_rows() %>%
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    #delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(iqr=paste0("(",iqr,")"),
         ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,strategy,quar_dur,`50%`,iqr,ui) %>% 
  htmlTable()

#ADHERENCE
plot_a_adherence <- get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    #adherence_iso==0.67,
    #adherence_quar==0.5,
    delay_scaling==1,
    !multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_iso,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               adherence_iso~adherence_quar, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence_quar =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 adherence_iso =
                   c("1" = "100% adhere\nto isolation",
                     "0.67" =
                       "67% adhere\nto isolation",
                     "0" =
                       "0% adhere\nto isolation")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_b_adherence <- get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    #adherence_iso==0.67,
    #adherence_quar==0.5,
    delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_iso,assay,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_iso,assay,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=strategy),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               adherence_iso~., labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence_quar =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 adherence_iso =
                   c("1" = "100% adhere\nto isolation",
                     "0.67" =
                       "67% adhere\nto isolation",
                     "0" =
                       "0% adhere\nto isolation")
               )) +
  plotting_theme

plot_a_adherence+plot_b_adherence+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "daily_vs_end_quar_adherence",
          base = "plot", 
          width = 300, 
          height = 200)


get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    #delay_scaling==0.5,
    !multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_iso,adherence_quar,assay,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_iso,adherence_quar,assay,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  mutate_at(.vars = vars(contains("%")), .funs = percent_format(accuracy = 1)) %>% 
  unite(iqr, c(`25%`,`75%`), sep = ", ") %>% 
  unite(ui, c(`2.5%`,`97.5%`), sep= ", ") %>% 
  mutate(iqr=paste0("(",iqr,")"),
         ui=paste0("(",ui,")")) %>% 
  arrange(-delay_scaling) %>% 
  select(delay_scaling,assay,strategy,quar_dur,`50%`,iqr,ui) %>% 
  htmlTable()

#by infection type
plot_a_type<- get(results_name) %>% 
  bind_rows() %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    !multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,adherence_quar,type,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,adherence_quar,assay,type,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               ~type, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_b_type <- get(results_name) %>% 
  bind_rows() %>% 
  filter(
    adherence_iso==0.67,
    adherence_quar==0.5,
    delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(ind_idx,strategy,assay,type,n_tests,delay_scaling,sampling_freq) %>% 
  filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>% 
  summarise(all=sum(inf_end-inf_start),
            prop=sum(max_overlap)/all
  ) %>% 
  group_by(strategy,assay,type,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Transmission potential averted")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               ~type, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme

plot_a_type+plot_b_type+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 400, 
          device = "png",
          prefix = "daily_vs_end_quar_type",
          base = "plot", 
          width = 300, 
          height = 150)
