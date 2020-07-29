source('packages.R')
source("utils.R")

main_scenarios <- input %>%  
  mutate(pre_board_screening = as.factor(pre_board_screening)) %>% 
  mutate(pre_board_screening = fct_explicit_na(pre_board_screening, "NA")) %>% 
  mutate(pre_board_screening = factor(pre_board_screening,
                                      levels = names(pre_board_labels),
                                      labels = pre_board_labels, ordered = T)) %>% 
  inner_join(main_scenarios) %>% 
  mutate(second_test_delay_ = 
           ifelse(is.na(second_test_delay),
                  0,
                  second_test_delay),
         time_in_iso = 
           first_test_delay + 
           second_test_delay_+
           post_flight_screening,
         time_in_iso = 
           first_test_delay+second_test_delay_+post_flight_screening) %>% 
  filter(!(!post_flight_screening & released_test == "Released after first test")) %>% 
  filter(!(post_flight_screening & released_test =="Released after mandatory isolation"))


arrival_released_person_days_summaries <- 
  arrival_released_times %>%
  filter(!is.na(days_released_inf), days_released_inf > 0) %>%
  complete(stage_released, released_test,sim, scenario, country) %>%
  mutate(days_released_inf = ifelse(is.na(days_released_inf), 0, days_released_inf)) %>%
  group_by(stage_released, released_test, scenario, country, sim) %>%
  summarise(days_released_inf = sum(days_released_inf)) %>%
  nest(data = -c(stage_released, released_test, scenario, country)) %>%
  mutate(Q = map(.x = data, ~quantile( .x$days_released_inf,
                                       probs = probs)),
         M = map_dbl(.x = data, ~mean(.x$days_released_inf))) %>%
  unnest_wider(Q) %>%
  select(-data)

figS5A_data <- input %>%
  mutate(second_test_delay_ = 
           ifelse(is.na(second_test_delay),
                  0,
                  second_test_delay),
         time_in_iso = 
           first_test_delay + 
           second_test_delay_+
           post_flight_screening) %>%
  
  inner_join(arrival_released_person_days_summaries) %>%
  mutate(pre_board_screening = as.factor(pre_board_screening)) %>% 
  mutate(pre_board_screening = fct_explicit_na(pre_board_screening, "NA")) %>% 
  mutate(pre_board_screening = factor(pre_board_screening,
                                      levels = names(pre_board_labels),
                                      labels = pre_board_labels, ordered = T)) %>% 
  inner_join(main_scenarios)  %>%  
  
  nest(data = -c(first_test_delay, second_test_delay)) %>%
  unite(col = "delays",
        first_test_delay, second_test_delay,
        sep = " + ", remove = FALSE) %>%
  unnest(data) %>%
  mutate(time_in_iso = factor(time_in_iso, 
                              levels = unique(.$time_in_iso),
                              ordered = T)) %>%
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
  filter(M!=0) %>%  
  mutate(released_test = factor(released_test,
                                levels = names(released_labels),
                                labels = released_labels, ordered = T)) 

figS5A <- figS5A_data %>% 
  filter(pre_board_screening == "None") %>%  
  make_release_figure(
    x = .,
    xlab = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab = ylabA,
    log_scale = FALSE,
    text_size = text_size,
    faceting = country ~ stringency)

arrival_released_person_days_summaries <- 
  arrival_released_times %>%
  filter(!is.na(days_released_inf), days_released_inf > 0) %>%
  complete(stage_released, released_test,sim, scenario, country) %>%
  mutate(days_released_inf = ifelse(is.na(days_released_inf), 0, days_released_inf)) %>%
  group_by(stage_released, released_test, scenario, country, sim) %>%
  summarise(days_released_inf = sum(days_released_inf)) 


baseline <- arrival_released_person_days_summaries %>% 
  filter(scenario==baseline_scenario) %>% 
  rename("baseline_days_released_inf"=days_released_inf,
         "baseline_scenario"=scenario,
         "baseline_released_test"=released_test) %>% 
  filter(stage_released=="Infectious",
         baseline_released_test=="Released after mandatory isolation") 

pd_risk_ratios <- arrival_released_person_days_summaries %>% 
  filter(stage_released=="Infectious",
         released_test=="Released after mandatory isolation"|
           released_test=="Released after first test"|
           released_test=="Released after second test") %>% 
  inner_join(baseline,by=c("stage_released", "sim", "country")) %>% 
  mutate(ratio=(days_released_inf)/(baseline_days_released_inf)) %>% 
  replace_na(list(ratio=1)) %>% 
  nest(data = -c(stage_released, released_test, scenario, country)) %>%
  mutate(Q = map(.x = data, ~quantile( .x$ratio,
                                       probs = probs)),
         M = map_dbl(.x = data, ~mean(.x$ratio))) %>%
  unnest_wider(Q) %>%
  select(-data)

figS5B_data <- main_scenarios %>% 
  inner_join(pd_risk_ratios) %>%
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

figS5B <- figS5B_data %>% 
  filter(pre_board_screening == "None") %>%  
  make_release_figure(
    x = .,
    xlab="Days in quarantine\n(including 1 day delay on testing results)",
    ylab="Ratio of cumulative infectious person-days remaining after release \nin comparison to low, no quarantine, no testing scenario",
    hline="dashed",
    text_size = text_size,
    log_scale=log_scale,
    faceting = country ~ stringency)+
  coord_cartesian(ylim=c(NA,1.2))

figS5 <- figS5A + figS5B + plot_layout(ncol = 1, guide = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))
