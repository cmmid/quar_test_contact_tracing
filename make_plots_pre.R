# source('packages.R')
# source("utils.R")
#covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

arrival_released_times_summaries <- 
  make_arrival_released_quantiles(arrival_released_times,
                                  country,
                                  pre_board_screening) 

figA_data <- plot_data(input, arrival_released_times_summaries, main_scenarios)

figS2A <- figA_data %>% 
  make_release_figure(
    x          = .,
    xlab       = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab       = ylabA,
    h_just     = 0,
    text_angle = 45,
    text_size  = 2.5,
    faceting   = country ~ stringency + pre_board_screening
    ) +
  theme(text = element_text(size = 18))

## DAYS OF INFECTIOUSNESS REMAINING

arrival_released_days_summaries <- 
  make_arrival_released_time_quantiles(arrival_released_times, 
                                       country,
                                       pre_board_screening)

figB_data <- plot_data(input, arrival_released_days_summaries, main_scenarios)

figS2B <- figB_data %>% 
  make_release_figure(
    x          = .,
    xlab       = "Days in quarantine\n(including 1 day delay on testing results)",
    ylab       = ylabB,
    h_just     =  0,
    text_angle = 45,
    text_size  =  2.5,
    faceting   = country ~ stringency + pre_board_screening) +
  theme(text = element_text(size = 18))+
  coord_cartesian(ylim = c(NA, 20))

figS2 <- figS2A + figS2B + plot_layout(ncol = 1, guide = "collect") +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))
