
# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("tracing_delays.R")
source("parameters.R")
source("kucirka_fitting.R")


seed <- 1000
n_ind_cases <- 1000
n_sec_cases <- 10
n_sims <- 100
asymp_parms <- asymp_fraction

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `10 day post-exposure quarantine only` = 
      data.frame(sampling_freq=NA,
                 tests=F,
                 multiple_tests=F,
                 assay=NA,
                 quar_dur=10),
    `Daily testing` = 
      data.frame(sampling_freq=1,
                 tests=T,
                 multiple_tests=T,
                 assay="LFA",
                 quar_dur=NA), 
    `10 day post-exposure quarantine with test` = 
      data.frame(sampling_freq=NA,
                 tests=T,
                 multiple_tests=F,
                 assay="LFA",
                 quar_dur = 10)) %>% 
      bind_rows(.id = "stringency")) %>% 
  crossing(post_symptom_window =  10,
           index_test_delay    =  c(1, 2, 3),  # time to entering quarantine
           delay_scaling       =  c(1, 0.5, 0),
           waning              =c("adhere_100"),
           test_sensitivity    =c(0.75,0.8,0.95),
           adherence=c(0,0.5,1)) %>% 
  mutate(test_to_tracing=3*delay_scaling) %>% 
  mutate(scenario=row_number()) %>% 
  #calculate time until release from exposure for each scenario
  mutate(time_since_exp=quar_dur)

# Filter input scenarios if required
input %<>% filter(
  index_test_delay %in% c(1),
  #delay_scaling    == 1,
  waning=="adhere_100",
  #quar_dur         %in% seq(0,14,by=2),
  #stringency       == "Test upon tracing\nand end of quarantine"
) %>% 
  mutate(scenario=row_number()) 

#browser()
set.seed(seed)



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
  crossing(distinct(input, index_test_delay, delay_scaling,test_to_tracing,test_sensitivity)) %>%     
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
                 delay_scaling,
                 test_sensitivity))

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


incubation_times_out <- left_join(input,
                                  incubation_times_out,
                                  by = c("index_test_delay", "delay_scaling","test_to_tracing")) %>% 
  mutate(adhering=rbinom(n=n(),size = 1,prob = adherence)) 

# generate testing times
incubation_times_out %<>% 
  mutate(test_t = pmap(.l = list(tracing_t=index_traced_t,
                                 sampling_freq=sampling_freq,
                                 max_time=10),
                      .f = test_times)) %>% 
  unnest(test_t)

#calc outcomes 
my_message("Calculating outcomes for each secondary case")
incubation_times_out %<>% calc_outcomes(x  = .,test_sensitivity=test_sensitivity)

#shift timings
incubation_times_out %<>% 
  mutate(exposed_q = 0 - 
           sec_onset_t + infect_shift,
         traced_q  = index_traced_t - 
           sec_onset_t + infect_shift,
         test_q    = test_t - 
           sec_onset_t + infect_shift, #isolate from point of test result
         onset_q = infect_shift) 

#find earliest positive test result
incubation_times_out %<>% 
  nest(test_t, test_p, test_q, test_no, test_label, screen) %>% 
  mutate(earliest_q      =  map_dbl(.f = earliest_pos, #for single tests - calc earliest positive test at which they isolate from
                                    .x = data)) %>% 
  select(-data)

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
    trans_pot_post_symp = case_when(
      type == "symptomatic" ~ pgamma(
        q = onset_q,
        shape = infect_shape,
        rate  = infect_rate,
        lower.tail = FALSE
      ),
      type == "asymptomatic" ~ 0
    )
  ) %>%
  mutate(
    trans_pot_pos_test = case_when(
      tests~pgamma(
      q = earliest_q,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
    ),
    TRUE~0
  )) %>%
  mutate(trans_pot_quar = case_when(
    !multiple_tests ~ pgamma(
      q = traced_q,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
    ),
    TRUE ~ 0
  )) %>%
  mutate(trans_pot_end_quar = case_when(
    !tests~pgamma(
      q= quar_dur - 
        sec_onset_t + infect_shift,
      shape = infect_shape,
      rate  = infect_rate,
      lower.tail = FALSE
    ),
    TRUE~0,
  )) %>% 
  mutate(trans_pot_quar=case_when(!tests~trans_pot_quar-trans_pot_end_quar,
                                  tests ~ trans_pot_quar-pmax(trans_pot_post_symp,trans_pot_pos_test))) %>% 
  mutate_at(
    .vars = vars(trans_pot_pos_test,
                 trans_pot_post_symp,
                 trans_pot_quar),
    .funs = function(x, y) {
      x / y
    },
    y = .$trans_pot_post_exp
  ) %>%
  mutate(
    trans_pot_post_symp =
      case_when(
        type == "symptomatic" & adhering == 1 ~
          trans_pot_post_symp,
        type == "symptomatic" & adhering == 0 ~
          0,
        TRUE ~ trans_pot_post_symp
      )
  ) %>%
  mutate(trans_pot_pos_test=case_when( 
         type == "symptomatic" &
           (trans_pot_pos_test > trans_pot_post_symp) ~
           trans_pot_pos_test - trans_pot_post_symp,
         type == "symptomatic" &
           (trans_pot_pos_test < trans_pot_post_symp) ~
           0,
         TRUE ~ trans_pot_pos_test
  )
) %>% 
  mutate(trans_pot_quar=case_when(
    #if you have a single test, you quarantine up until your positive test or symptom onset; if positive, you self-isolate longer
    !multiple_tests&tests&trans_pot_quar>(trans_pot_post_symp|trans_pot_pos_test)~
      trans_pot_quar-pmax(trans_pot_post_symp,trans_pot_pos_test),
    !multiple_tests&tests&trans_pot_quar<trans_pot_post_symp~
      0,
    #if you don't have any tests, you quarantine until 10 days or symptom onset
    !multiple_tests&!tests&trans_pot_quar>trans_pot_post_symp~
      trans_pot_quar-trans_pot_post_symp,
    #if you have multiple tests
    TRUE ~ 0
  )) %>% 
  mutate(
    trans_pot_averted = trans_pot_quar+trans_pot_pos_test+trans_pot_post_symp) 


incubation_times_out %>% left_join(input) %>% 
  summarise_simulation(faceting=stringency~type+delay_scaling+adherence) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate(stringency=fct_relevel(stringency, "10 day post-exposure quarantine only", "10 day post-exposure quarantine with test", "Daily testing")) %>% 
  filter(yvar=="trans_pot_averted",
         type=="all",
         adherence==0.5) %>% 
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
  scale_x_discrete(labels=delay_scaling_labeller)+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  # )+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1),breaks = breaks_width(0.25))+
  labs(x="",
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
               adherence ~ ., labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% self-isolate\nupon symptom onset",
                     "0.5" =
                       "50% self-isolate\nupon symptom onset",
                     "0" =
                       "0% self-isolate\nupon symptom onset")
               )) + 
  plotting_theme+
  scale_colour_brewer(name="Strategy",palette = "Dark2")+
  scale_fill_brewer(name="Strategy",palette="Dark2")

save_plot(dpi = 400, 
                    device = "png",
                    prefix = "daily_vs_end_quar",
                    base = "plot", 
                    #width = 500, 
                    height = 125)

incubation_times_out %>% left_join(input) %>% 
  summarise_simulation(faceting=stringency~type+delay_scaling+adherence+test_sensitivity) %>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate(stringency=fct_relevel(stringency, "10 day post-exposure quarantine only", "10 day post-exposure quarantine with test", "Daily testing")) %>% 
  filter(yvar=="trans_pot_averted",
         type=="all",
        # adherence==0.5
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
               test_sensitivity ~ adherence, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% self-isolate\nupon symptom onset",
                     "0.5" =
                       "50% self-isolate\nupon symptom onset",
                     "0" =
                       "0% self-isolate\nupon symptom onset")
               )) + 
  plotting_theme+
  scale_colour_brewer(name="Strategy",palette = "Dark2")+
  scale_fill_brewer(name="Strategy",palette="Dark2")

save_plot(dpi = 400, 
                    device = "png",
                    prefix = "daily_vs_end_quar_all",
                    base = "plot", 
                    #width = 500, 
                    height = 150)
