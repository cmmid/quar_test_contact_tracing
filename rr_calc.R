#source('packages.R')
#source("utils.R")

# THIS WHOLE THING NEEDS A REWRITE
# Should be set up as a function
# pass in the output of a simulation
# pass in a data.frame containing details of the scenario to use as baseline
## error checking to ensure that joining on the input baseline doesn't result
## in a data fram that is "too long"

main_scenarios <- input %>% 
  inner_join(main_scenarios) %>% 
  mutate(second_test_delay_ = 
           ifelse(is.na(second_test_delay),
                  0,
                  second_test_delay),
         time_in_iso = 
           first_test_delay + 
           second_test_delay_+
           screening,
         time_in_iso = 
           first_test_delay+second_test_delay_+screening) %>% 
  filter(!(!screening & released_test == "Released after first test")) %>% 
  filter(!(screening & released_test =="Released after mandatory isolation"))

released_times_summaries <- 
  results %>% 
  mutate(time_in_iso=released_t-traced_t) %>% 
  count(stage_released, released_test, sim, scenario) %>%
  complete(stage_released, released_test, sim, scenario) %>% 
  mutate(n=replace_na(n,0))


baseline <- released_times_summaries %>% 
  filter(scenario==21) %>% 
  rename("baseline_n"=n,
         "baseline_scenario"=scenario,
         "baseline_released_test"=released_test) %>% 
  filter(stage_released=="Infectious",
         baseline_released_test=="Released after mandatory isolation") 

n_risk_ratios <- released_times_summaries %>% 
  filter(stage_released=="Infectious",
         released_test=="Released after mandatory isolation"|
           released_test=="Released after first test"|
           released_test=="Released after second test") %>% 
  inner_join(baseline) %>% 
  mutate(ratio=(n)/(baseline_n)) %>% 
  replace_na(list(ratio=1)) %>% 
  nest(data = -c(stage_released, released_test, scenario,baseline_scenario)) %>%
  mutate(Q = map(.x = data, ~quantile(.x$ratio, probs = probs)),
         M = map_dbl(.x = data, ~mean(.x$ratio))) %>%
  unnest_wider(Q) %>%
  select(-data)

n_fig_data <- main_scenarios %>% 
  inner_join(n_risk_ratios) %>%
  nest(data = -c(first_test_delay, second_test_delay)) %>%
  unite(col = "delays",
        first_test_delay, second_test_delay,
        sep = " + ", remove = FALSE) %>%
  unnest(data) %>%
  mutate(time_in_iso = factor(time_in_iso, 
                              levels = unique(.$time_in_iso),
                              ordered = T)) %>%
  mutate(time_in_iso=fct_reorder(time_in_iso,as.numeric(as.character(time_in_iso)))) %>% 
  mutate(stringency = factor(stringency,
                             levels = c("low",
                                        "moderate",
                                        "high",
                                        "maximum"),
                             labels = c("Low",
                                        "Mod.",
                                        "High",
                                        "Max."),
                             ordered = T)) %>%
  mutate(released_test = factor(released_test,
                                levels = names(released_labels),
                                labels = released_labels, ordered = T)) 


ggplot(n_fig_data)+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_pointrange(aes(x=time_in_iso,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=released_test,group=delays),position = position_dodge2(width=0.75))+
  scale_y_log10()+
  facet_grid(~stringency,scales="free",space="free")+
  scale_colour_manual(values=covid_pal)+
  theme_minimal()+
  theme(panel.border = element_rect(fill = NA))

# n_fig <- n_fig_data %>% 
#   make_release_figure(
#     x = .,
#     xlab="Days in quarantine\n(including 1 day delay on testing results)",
#     ylab=ifelse(baseline_scenario=="87","Ratio of infectious persons released per week \nin comparison to high, 14 day quarantine, no testing scenario","Ratio of per-individual infectious persons released per week \nin comparison to low, no quarantine, no testing scenario"),
#     hline = "dashed",
#     text_size = text_size,
#     log_scale = log_scale,
#     faceting = ~ stringency)

released_person_days_summaries <- 
  results %>%
  filter(!is.na(days_released_inf),days_released_inf > 0) %>%
  complete(stage_released, released_test,sim, scenario) %>%
  mutate(days_released_inf=replace_na(days_released_inf,0)) 

baseline <- released_person_days_summaries %>% 
  filter(scenario==21) %>% 
  rename("baseline_days_released_inf"=days_released_inf,
         "baseline_scenario"=scenario,
         "baseline_released_test"=released_test) %>% 
  filter(stage_released=="Infectious",
         baseline_released_test=="Released after mandatory isolation") %>% 
  select(c(stage_released,baseline_released_test,sim,idx,baseline_scenario,baseline_days_released_inf))

pd_risk_ratios <- released_person_days_summaries %>% 
  filter(stage_released=="Infectious",
         released_test=="Released after mandatory isolation"|
           released_test=="Released after first test"|
           released_test=="Released after second test") %>% 
  inner_join(baseline,by=c("stage_released","sim","idx")) %>% 
  mutate(ratio=(days_released_inf)/(baseline_days_released_inf)) %>% 
  replace_na(list(ratio=1)) %>% 
  nest(data = -c(stage_released, released_test, scenario)) %>%
  mutate(Q = map(.x = data, ~quantile( .x$ratio,
                                       probs = probs)),
         M = map_dbl(.x = data, ~mean(.x$ratio))) %>%
  unnest_wider(Q) %>%
  select(-data)

pd_fig_data <-pd_risk_ratios  %>% 
  inner_join(main_scenarios) %>%
  nest(data = -c(first_test_delay, second_test_delay)) %>%
  unite(col = "delays",
        first_test_delay, second_test_delay,
        sep = " + ", remove = FALSE) %>%
  unnest(data) %>%
  mutate(time_in_iso = factor(time_in_iso, 
                              levels = unique(.$time_in_iso),
                              ordered = T)) %>%
  mutate(time_in_iso=fct_reorder(time_in_iso,as.numeric(as.character(time_in_iso)))) %>% 
  mutate(stringency = factor(stringency,
                             levels = c("low",
                                        "moderate",
                                        "high",
                                        "maximum"),
                             labels = c("Low",
                                        "Mod.",
                                        "High",
                                        "Max."),
                             ordered = T)) %>%
  mutate(released_test = factor(released_test,
                                levels = names(released_labels),
                                labels = released_labels, ordered = T)) 


ggplot(pd_fig_data)+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_pointrange(aes(x=time_in_iso,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=released_test,group=delays),position = position_dodge2(width=0.75))+
  scale_y_log10()+
  facet_grid(~stringency,scales="free",space="free")+
  scale_colour_manual(values=covid_pal)+
  theme_minimal()+
  theme(panel.border = element_rect(fill = NA))


rr_figs <- n_fig + pd_fig + plot_layout(ncol = 1, guide = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))
