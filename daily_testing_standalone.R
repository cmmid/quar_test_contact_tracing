
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("kucirka_fitting.R")
source("parameters.R")


input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
     `Daily testing` = 
      crossing(sampling_freq=1,
                 tests=T,
                 multiple_tests=T,
                 assay="LFA",
                 n_tests = c(1,3,5,7,10), 
                 quar_dur=NA), 
     `Post-exposure quarantine only` = 
      crossing(sampling_freq=NA,
                 tests=F,
                 multiple_tests=F,
                 assay=NA,
                 n_tests=NA,
                 quar_dur=c(0,3,5,7,10,14)),
    `Post-exposure quarantine with test` = 
      crossing(sampling_freq=NA,
                 tests=T,
                 multiple_tests=F,
                 n_tests=NA,
                 assay="LFA",
                 quar_dur=c(0,3,5,7,10,14))) %>% 
      bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window =  10,
           index_test_delay    =  c(1),  # time to entering quarantine (index cases)
           delay_scaling       =  c(1, 0.5),
           test_sensitivity    =c(0.75,0.8,0.9),
           adherence_quar=c(0,0.5,1),
           adherence_iso=c(0,0.67,1)) %>% 
  mutate(test_to_tracing=3*delay_scaling) %>% 
  #filter(tests&!multiple_tests) %>% 
  mutate(scenario=row_number()) 

input_split <-
  input %>% 
  rowwise %>%
  group_split

run_model <- function(
  n_sims          = 1000,
  n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
  n_ind_cases     = 10000,
  input,
  seed            = 145,
  asymp_parms,
  return_full = TRUE,
  faceting = stringency ~ type){
  
 #browser()
  
  set.seed(seed)
  
 message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), input$scenario))
  
#Generate incubation periods to sample
incubation_times <- make_incubation_times(
  n_travellers = n_ind_cases,
  pathogen     = pathogen,
  asymp_parms  = asymp_parms)



inf <- data.frame(prop_asy = rbeta(n = n_sims,
                                   shape1 = asymp_parms$shape1,
                                   shape2 = asymp_parms$shape2)) 



# Generate index cases' inc times
ind_inc <- incubation_times %>% 
  filter(type=="symptomatic") %>% 
  sample_n(n_sims) %>% 
  mutate(sim = seq(1L, n_sims, by = 1L)) %>% 
  bind_cols(inf) %>% 
  #sample test result delay
  ## sample uniformly between 0 and 1 when 0.5...
  crossing(distinct(input, index_test_delay, delay_scaling,test_to_tracing,test_sensitivity)) %>%     
  rename("index_onset_t" = onset) %>% 
  mutate(index_testing_t    = index_onset_t + index_test_delay,
         sec_traced_t     = index_onset_t + index_test_delay + test_to_tracing)

#rm(list = c("P_t", "P_r", "P_c", "inf"))

my_message("Generating secondary cases' incubation times")

#Generate secondary cases
sec_cases <- make_incubation_times(
  n_travellers = n_sec_cases,
  pathogen     = pathogen,
  asymp_parms  = asymp_parms)

my_message("Generating secondary cases' exposure times")
ind_inc %<>% 
  nest(data = -c(sim,
                 prop_asy,
                 test_to_tracing,
                 index_onset_t,
                 index_test_delay,
                 index_testing_t,
                 sec_traced_t,
                 delay_scaling,
                 test_sensitivity))

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

rm(ind_inc)

incubation_times_out <- left_join(input,
                                  incubation_times_out,
                                  by = c("index_test_delay", "delay_scaling","test_to_tracing","test_sensitivity")) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso)) 

#browser()

# generate testing times
incubation_times_out %<>% 
  mutate(test_t = pmap(.l = list(tracing_t=sec_traced_t,
                                 sampling_freq=sampling_freq,
                                 max_time=quar_dur,
                                 max_tests=n_tests),
                      .f = test_times)) %>% 
  unnest(test_t)

#calc outcomes 
my_message("Calculating outcomes for each secondary case")
incubation_times_out %<>% calc_outcomes(x  = .,test_sensitivity=test_sensitivity)

#shift timings
incubation_times_out %<>% 
  mutate(exposed_q = sec_exposed_t - 
           sec_onset_t + infect_shift,
         traced_q  = sec_traced_t - 
           sec_onset_t + infect_shift,
         test_q    = test_t - 
           sec_onset_t + infect_shift, #isolate from point of test result
         onset_q = infect_shift,
         quar_end_q= exposed_q + quar_dur)

#browser()

#find earliest positive test result
incubation_times_out %<>% 
  nest(test_t, test_p, test_q, test_no, test_label, screen) %>% 
  mutate(earliest_q      =  map(.f = earliest_pos2, 
                                    .x = data)) %>% 
  unnest_wider(earliest_q) %>% 
  rename("earliest_q"=test_q)

# calculate remaining transmission potential averted by positive test
incubation_times_out %<>%
  mutate(
    trans_pot_post_exp =
      pgamma(
        q     = exposed_q,
        shape = infect_shape,
        rate  = infect_rate,
        lower.tail = FALSE
      )
  ) %>%
  mutate(
    trans_pot_start_symp = case_when(
      type == "symptomatic"  ~ pgamma(
        q = onset_q,
        shape = infect_shape,
        rate  = infect_rate,
        lower.tail = FALSE
      )*adhering_iso,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    trans_pot_end_symp = case_when(
      type == "symptomatic"  ~ pgamma(
        q = onset_q+10,
        shape = infect_shape,
        rate  = infect_rate,
        lower.tail = FALSE
      )*adhering_iso,
      TRUE ~ 0
    )
  ) %>%
  mutate(trans_pot_symp=(trans_pot_start_symp-trans_pot_end_symp)) %>% 
  mutate(
    trans_pot_pos_test = case_when(
      tests~pgamma(
      q = earliest_q,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
    )*adhering_iso,
    TRUE~0
  )) %>%
  mutate(trans_pot_end_test = case_when(
    tests~ pgamma(
      q = earliest_q+10,
     shape = infect_shape,
     rate  = infect_rate,
     lower.tail = FALSE
  )*adhering_iso,
  TRUE~0)) %>% 
  mutate(trans_pot_test = trans_pot_pos_test-trans_pot_end_test) %>% 
  mutate(trans_pot_traced =  pgamma(
      q = traced_q,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
  )*adhering_quar) %>%
  mutate(trans_pot_end_quar = case_when(
    !multiple_tests~pgamma(
      q= quar_end_q,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
    )*adhering_quar,
    TRUE~0,
  )) %>% 
  mutate(trans_pot_quar=case_when(!tests~(trans_pot_traced-trans_pot_end_quar),
                                  tests&!multiple_tests~(trans_pot_traced-trans_pot_end_quar)+trans_pot_test,
                                  multiple_tests~0)) %>% 
  mutate_at(
    .vars = vars(trans_pot_pos_test,
                 trans_pot_end_test,
                 trans_pot_test,
                 trans_pot_start_symp,
                 trans_pot_symp,
                 trans_pot_traced,
                 trans_pot_end_quar,
                 trans_pot_quar),
    .funs = function(x, y) {
      x / y
    },
    y = .$trans_pot_post_exp
  ) %>%
  mutate(trans_pot_averted=pmax(trans_pot_symp,trans_pot_test,trans_pot_quar))
#browser()
}

results_name <- "results_list"

assign(x     = results_name,
       value = map(
         .x =  input_split,
         .f =  ~ run_model(
           input=.x,
           seed = 1000,
           n_ind_cases = 1000,
           n_sec_cases = 10,
           n_sims = 100,
           asymp_parms = asymp_fraction
         )))




 saveRDS(get(results_name),"daily_results.rds")
# incubation_times_out <- readRDS("daily_results.rds")




get(results_name) %>% 
  bind_rows() %>% 
  filter(test_sensitivity==0.75

    ) %>%
  group_by(sim,stringency,type,adherence_quar,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  summarise(n=n(),
            prop=sum(trans_pot_averted)/n) %>% 
  group_by(stringency,type,adherence_quar,quar_dur,n_tests,test_sensitivity,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency=factor(stringency)) %>% 
ggplot(aes(x = factor(delay_scaling), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.2),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.2),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.2)) +
  scale_x_discrete(labels=delay_scaling_labeller,guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  # )+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x="",
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               type ~ test_sensitivity+quar_dur+n_tests, labeller = labeller(
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
  scale_colour_brewer(name="Strategy",palette = "Dark2")+
  scale_fill_brewer(name="Strategy",palette="Dark2")

save_plot(dpi = 400, 
                    device = "png",
                    prefix = "daily_vs_end_quar_n_tests",
                    base = "plot", 
                    #width = 500, 
                    height = 125)

incubation_times_out %>% left_join(input) %>% 
  summarise_simulation(faceting=stringency~type+delay_scaling+adherence+test_sensitivity) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate(stringency=fct_relevel(stringency, "10 day post-exposure quarantine only", "10 day post-exposure quarantine with test", "Daily testing")) %>% 
  filter(yvar=="trans_pot_averted",
         type=="all",
         test_sensitivity==0.75,
         #adherence==0.5
        ) %>% 
  ggplot(aes(x = factor(delay_scaling), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.2,
             position=position_dodge(width=0.5)) +
  scale_x_discrete(labels=delay_scaling_labeller,guide=guide_axis(angle = 45))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  # )+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x="",
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
                ~ adherence, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 test_sensitivity=c("0.75"="0.75 * PCR curve",
                                    "0.8"="0.8 * PCR curve",
                                    "0.9"="0.9 * PCR curve")
                 
               )) + 
  plotting_theme+
  scale_colour_brewer(name="Strategy",palette = "Dark2")+
  scale_fill_brewer(name="Strategy",palette="Dark2")

save_plot(dpi = 400, 
                    device = "png",
                    prefix = "daily_vs_end_quar_all",
                    base = "plot", 
                    width = 250, 
                    height = 150)

incubation_times_out %>% left_join(input) %>% 
  summarise_simulation(faceting=stringency~type+delay_scaling+adherence+test_sensitivity) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate(stringency=fct_relevel(stringency, "10 day post-exposure quarantine only", "10 day post-exposure quarantine with test", "Daily testing")) %>% 
  filter(yvar=="trans_pot_averted",
         type=="all",
         #test_sensitivity==0.75,
         adherence==0.5
  ) %>% 
  ggplot(aes(x = factor(delay_scaling), y = `50%`)) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=stringency),position=position_dodge(width=0.5),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=stringency),
             #pch="-",
             size=1.2,
             position=position_dodge(width=0.5)) +
  scale_x_discrete(labels=delay_scaling_labeller,guide=guide_axis(angle = 45))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  # )+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x="",
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               ~ adherence+test_sensitivity, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine"),
                 test_sensitivity=c("0.75"="0.75 * PCR curve",
                                    "0.8"="0.8 * PCR curve",
                                    "0.9"="0.9 * PCR curve")
                 
               )) + 
  plotting_theme+
  scale_colour_brewer(name="Strategy",palette = "Dark2")+
  scale_fill_brewer(name="Strategy",palette="Dark2")

save_plot(dpi = 400, 
          device = "png",
          prefix = "daily_vs_end_quar_test_sens",
          base = "plot", 
          width = 250, 
          height = 150)
