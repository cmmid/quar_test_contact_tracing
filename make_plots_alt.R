# source('packages.R')
# source("utils.R")


arrival_released_times_summaries <- 
  make_arrival_released_quantiles(arrival_released_times, country)


fig3A_data <- input %>%
  mutate(second_test_delay_ = 
           ifelse(is.na(second_test_delay),
                  0,
                  second_test_delay),
         time_in_iso = 
           first_test_delay + 
           second_test_delay_+
           post_flight_screening) %>% 
  inner_join(arrival_released_times_summaries) %>%
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


fig3A <- fig3A_data %>% 
  filter(pre_board_screening == "None") %>% 
  make_release_figure(
    x = .,
    xlab = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab = ylabA,
    text_size = text_size,
    faceting = country ~ stringency)

## person-days


arrival_released_person_days_summaries <- 
  arrival_released_times %>%
  filter(!is.na(days_released_inf),days_released_inf>0) %>%
  complete(stage_released, released_test,sim,scenario, country) %>%
  mutate(days_released_inf=replace_na(days_released_inf,0)) %>% 
  #filter(days_released_inf > 0) %>% 
  nest(data = -c(stage_released, released_test, scenario, country)) %>%
  mutate(Q = map(.x = data, ~quantile( .x$days_released_inf,
                                       probs = probs)),
         M = map_dbl(.x = data, ~mean(.x$days_released_inf))) %>%
  unnest_wider(Q) %>%
  select(-data)

fig3B_data <- input %>%
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

fig3B <- fig3B_data %>% 
  filter(pre_board_screening == "None") %>%  
  make_release_figure(
    x = .,
    xlab = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab = ylabB,
    text_size = text_size,
    faceting = country ~ stringency)+
  coord_cartesian(ylim=c(NA,16))

fig3 <- fig3A + fig3B + plot_layout(ncol = 1, guide = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))
