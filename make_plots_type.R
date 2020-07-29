#source('packages.R')
#source("utils.R")

arrival_released_times_summaries <- 
  make_arrival_released_quantiles(arrival_released_times, country, type)

figS3A_data <- plot_data(input, arrival_released_times_summaries, main_scenarios)

figS3A <- figS3A_data %>% 
  filter(pre_board_screening == "None") %>% 
  make_release_figure(
    x         = .,
    text_size = text_size,
    xlab      = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab      = ylabA) +
  facet_nested(country + type ~ stringency, 
               nest_line = T,
               labeller  = labeller(type=type_labels),
               scales    = "free_x",
               space     = "free_x")

## person-days


arrival_released_days_summaries <- 
  make_arrival_released_time_quantiles(arrival_released_times, country, type)

figS3B_data <- plot_data(input, arrival_released_days_summaries, main_scenarios)

figS3B <- figS3B_data %>% 
  filter(pre_board_screening == "None") %>%  
  make_release_figure(
    x = .,
    xlab = "Days in quarantine\n(including 1 day delay on testing results)",
    text_size= text_size,
    ylab = ylabB)+
    facet_nested(country+type ~ stringency, 
                 nest_line = T,
               scales = "free_x",space = "free_x",
               labeller = labeller(type=type_labels))+
  coord_cartesian(ylim=c(NA,20))


figS3 <- figS3A + figS3B + plot_layout(ncol = 1, guide = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))
