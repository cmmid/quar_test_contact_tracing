source('packages.R')
source("utils.R")


arrival_released_times <- 
  #readRDS("data/arrival_released_times_200_sims_20200713_175927.RDS") %>%
  arrival_released_times %>% 
  mutate(pre_board_screening = as.factor(pre_board_screening)) %>% 
  mutate(pre_board_screening = fct_explicit_na(pre_board_screening, "NA")) %>% 
  mutate(pre_board_screening = factor(pre_board_screening,
                                      levels = names(pre_board_labels),
                                      labels = pre_board_labels, ordered = T))


source('kucirka_fitting.R')
kucirka_plot <- ggplot(data = dat_pred, aes(x = day, y = fit)) + 
  # geom_ribbon(aes(ymin = L,
  #                 ymax = U),
  #             fill = lshtm_greens[2],
  #             alpha = 0.25) + 
  geom_line(color = lshtm_greens[2]) +
  geom_point(data = dat, aes(y =pct_pos),
             pch = 16,
             color = lshtm_greens[1],
             position = position_jitter(w = 0.05, h=0.005),
             alpha = 0.25) +
  theme_bw() +
  xlab("Days past presumed exposure") +
  ylab("Proportion of SARS-CoV-2\ninfections detected")


distributions_plot <- 
  incubation_times %>%
  transmute(`Time to onset of symptoms`      = exp_to_onset,
            `End of symptomatic period`      = onset_to_recov + exp_to_onset,
            `Duration of symptomatic period` = symp_dur,
            `Start of infectious period`     = inf_start,
            `End of infectious period`       = inf_end,
            `Duration of infectious period`  = inf_dur,
            type = type) %>%
  gather(key, value, -type) %>%
  mutate(key = fct_inorder(key),
         type = str_to_title(type),
         type = fct_rev(factor(type))) %>%
  filter(!(type == "Asymptomatic" & grepl(pattern = "symptom", x = key))) %>%
  mutate(value  = pmin(value, 31)) %>%
  ggplot(data = .,
         aes(x = value)) +
  geom_histogram(binwidth = 1, center = 0.5,
                 aes(y = ..density..),
                 color = "black", fill = lshtm_greens[1]) +
  facet_wrap(~ type + key, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 30, by = 10),
                     limits = c(0, 31),
                     labels = c("0", "10", "20", "30+")) +
  xlab("Time (days)") + 
  ylab("Density") +
  theme_minimal() +
  theme(panel.border = element_rect(fill=NA),
        #text = element_text(size=20),
        panel.spacing.x = unit(12,"mm"),
        panel.spacing.y = unit(6,"mm"))


fig1 <- patchwork::wrap_plots(kucirka_plot, 
                              distributions_plot, ncol = 1, heights = c(1,2))

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig1_.",.x),
              plot     = fig1,
              width    = 210,
              height   = 210,
              units    = "mm",
              dpi      = 300))
### end distributions


arrival_released_times_summaries <- 
  make_arrival_released_quantiles(arrival_released_times, country) 

plot_data <- input %>%
  mutate(second_test_delay_ = 
           ifelse(is.na(second_test_delay),
                  0,
                  second_test_delay),
         time_in_iso = 
           first_test_delay + 
           second_test_delay_+
           post_flight_screening) %>% 
  mutate(pre_board_screening = as.factor(pre_board_screening)) %>% 
  mutate(pre_board_screening = fct_explicit_na(pre_board_screening, "NA")) %>% 
  mutate(pre_board_screening = factor(pre_board_screening,
                                      levels = names(pre_board_labels),
                                      labels = pre_board_labels, ordered = T)) %>% 
  inner_join(arrival_released_times_summaries) %>%
  inner_join(main_scenarios)  %>%
  nest(data = -c(first_test_delay, second_test_delay)) %>%
  unite(col = "delays",
        first_test_delay, second_test_delay,
        sep = " + ", remove = FALSE) %>%
  unnest(data) %>%
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
  mutate(released_test = fct_rev(fct_recode(released_test,
                                    "No post-flight test"  = "Released after mandatory isolation",
                                    "One post-flight test" = "Released after first test"))) %>% 
   mutate(time_in_iso = factor(time_in_iso, 
                                  levels = unique(.$time_in_iso),
                                  ordered = T)) %>% 
  mutate(time_in_iso=as.numeric(as.character(time_in_iso))) %>% 
  filter(M!=0)

### Days in quarantine vs n released ----

main_plot <- ggplot()+
  geom_linerange(data=plot_data %>% 
                   filter(stringency!="High") %>% 
                   filter(pre_board_screening=="No pre-flight test"),
                 aes(x=time_in_iso,
                     ymin=`2.5%`,
                     ymax=`97.5%`,
                     colour=released_test),
                 size=3,
                 alpha=0.6)+
  geom_linerange(data=plot_data %>% 
                   filter(stringency!="High")%>% 
                   filter(pre_board_screening=="No pre-flight test"),
                 aes(x=time_in_iso,
                     ymin=`25%`,
                     ymax=`75%`,
                     colour=released_test),
                 size=3)+
  scale_colour_manual(name="",
                      values=lshtm_greens)+
  ggnewscale::new_scale_colour()+
  geom_linerange(data=plot_data %>% 
                   filter(stringency=="High")%>% 
                   filter(pre_board_screening=="No pre-flight test")%>% 
                   left_join(bivariate_color_scale) ,
                 aes(x=time_in_iso,
                     ymin=`2.5%`,
                     ymax=`97.5%`,
                     colour=colour),
                 position = position_dodge2(width=0.6),
                 size=3,
                 alpha=0.6)+
  geom_linerange(data=plot_data %>% 
                   filter(stringency=="High")%>% 
                   filter(pre_board_screening=="No pre-flight test") %>% 
                   left_join(bivariate_color_scale) ,
                 aes(x=time_in_iso,
                     ymin=`25%`,
                     ymax=`75%`,
                     colour=colour),
                 position = position_dodge2(width=0.6),
                 size=3)+
  geom_point(data=plot_data %>%
               filter(stringency!="High") %>%
               filter(pre_board_screening=="No pre-flight test") ,
             aes(x=time_in_iso,y = M, group = released_test),
             color="white",
             shape=1,
             fill=NA,
             size=2)+
  geom_point(data=plot_data %>%
               filter(stringency=="High") %>%
               filter(pre_board_screening=="No pre-flight test")%>%
               left_join(bivariate_color_scale),
             aes(x=time_in_iso,y = M, group = colour),
             color="white",
             shape=1,
             fill=NA,
             size=2,
             position = position_dodge2(width = 0.6))+
  scale_shape(solid=FALSE)+
  scale_colour_identity(guide=FALSE)+
  scale_x_continuous(breaks = breaks_width(1))+
  scale_y_continuous(limits=c(0,NA),
                     breaks = breaks_width(5))+
  coord_cartesian(expand = TRUE)+
  facet_grid(country~stringency,
               scales="free_x",
               space = "free_x",
               labeller=labeller(stringency=capitalize)) +
  #force_panelsizes(rows=8,cols=c(1.5,6,8,1.5))+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        text = element_text(size=20),
        # panel.spacing.x = unit(12,"mm"),
        # panel.spacing.y = unit(6,"mm"),
        legend.position = c(0.35,0.95),
        strip.placement = "outside")+
  xlab("Days in quarantine")+
  ylab("Number of infectious persons released per week")

fig2 <- ggdraw()+
  draw_plot(main_plot,0,0,1,1)+
  draw_plot(legend,0.6,0.73,0.3,0.2)
fig2

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig2.",.x),
              plot=fig2,
              width = 420, height = 280, units="mm",
              dpi = 320))
# Person days spent infectious after release ----
days_plot_data <- arrival_released_times %>%
  mutate(time_in_iso=released_t-flight_arrival) %>%
  inner_join(main_scenarios) %>%
  unite("delays",first_test_delay,second_test_delay,sep = " + ") %>%
  mutate_at(.vars=vars(sim,
                       stringency,
                       country,
                       time_in_iso,
                       delays,
                       pre_board_screening,
                       post_flight_screening,
                       released_test),
            .funs = as.factor) %>%
  group_by(sim,
           stringency,
           country,
           time_in_iso,
           delays,
           pre_board_screening,
           post_flight_screening,
           released_test,.drop = F) %>%
  summarise(sum_days_inf=sum(days_released_inf)) %>%
  group_by_at(.vars = vars(stringency,
                           country,
                           time_in_iso,
                           delays,
                           pre_board_screening,
                           post_flight_screening,
                           released_test)) %>%
  summarise(qs = quantile(sum_days_inf, probs),
            probs = probs) %>%
  pivot_wider(names_from = probs,values_from=qs) %>%
  mutate(sum=`0.025`+`0.25`+`0.5`+`0.75`+`0.975`) %>%
  filter(sum!=0) %>%
  mutate(test_facet = case_when(stringency == "low"|
                                  stringency =="moderate" |
                                  stringency =="maximum" ~ "No tests or one test",
                                stringency == "high" ~ "Two tests")) %>% 
  mutate(time_in_iso=as.numeric(as.character(time_in_iso))) %>% 
  mutate(released_test = fct_rev(fct_recode(released_test,
                                    "No test"  = "Released after mandatory isolation",
                                    "One test" = "Released after first test"))) %>% 
  mutate(stringency=fct_relevel(as.factor(stringency),"low","moderate","high","maximum",after=Inf)) 

days_plot <- ggplot()+
  geom_linerange(data=days_plot_data %>% 
                   filter(stringency!="high") %>% 
                   filter(pre_board_screening=="No pre-flight test") ,
                 aes(x=time_in_iso,
                     ymin=`0.025`,ymax=`0.975`,
                     colour=released_test),
                 position = position_dodge2(width=0.6),
                 size=2,
                 alpha=0.7)+
  geom_linerange(data=days_plot_data %>% 
                   filter(stringency!="high") %>% 
                   filter(pre_board_screening=="No pre-flight test"),
                 aes(x=time_in_iso,
                     ymin=`0.25`,
                     ymax=`0.75`,
                     colour=released_test),
                 position = position_dodge2(width=0.6),
                 size=2)+
  scale_colour_manual(name="",
                      values=lshtm_greens,guide=F)+
  ggnewscale::new_scale_colour()+
  geom_linerange(data=days_plot_data %>% 
                   filter(stringency=="high") %>% 
                   filter(pre_board_screening=="No pre-flight test") %>% 
                   left_join(bivariate_color_scale),
                 aes(x=time_in_iso,
                     ymin=`0.025`,
                     ymax=`0.975`,
                     colour=colour),
                 position = position_dodge2(width=0.6),
                 size=2,
                 alpha=0.7)+
  geom_linerange(data=days_plot_data %>% 
                   filter(stringency=="high") %>% 
                   filter(pre_board_screening=="No pre-flight test") %>% 
                   left_join(bivariate_color_scale),
                 aes(x=time_in_iso,
                     ymin=`0.25`,
                     ymax=`0.75`,
                     colour=colour),
                 position = position_dodge2(width=0.6),
                 size=2)+
  scale_colour_identity(guide=FALSE)+
  scale_x_continuous(breaks=breaks_width(1))+
  scale_y_continuous(limits=c(0,NA),
                     #breaks = breaks_width(20),
                     expand = expansion(mult = c(0, .1)))+
  facet_nested(country~stringency,scales="free_x",
               labeller=labeller(stringency=capitalize)) + 
  force_panelsizes(rows=8,cols=c(2,6,8,2))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA),
        panel.spacing.x = unit(12,"mm"),
        panel.spacing.y = unit(6,"mm"),
        text = element_text(size=20),
        legend.position = c(0.35,0.95),
        strip.placement = "outside")+
  xlab("Days in quarantine")+
  ylab("Cumulative infectious person-days remaining after release")

fig4 <- ggdraw()+
  draw_plot(days_plot,0,0,1,1)


list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig4.",.x),
              plot=fig4,
              width = 420, height = 280, units="mm",
              dpi = 320))

fig2ab <- cowplot::plot_grid(fig2,NULL,fig4, rel_heights=c(1,0.1,1),ncol=1,labels = c("A","","B"))

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig2ab.",.x),
              plot=fig2ab,
              width = 420, height = 550, units="mm",
              dpi = 200))

# Pre-board screening ----

pre_board_plot <- ggplot()+
  ggnewscale::new_scale_colour()+
  geom_linerange(data=plot_data %>% mutate(stringency=fct_recode(stringency,"max"="maximum")) %>% 
                   mutate(stringency=fct_relevel(as.factor(stringency),"low","moderate","high","max",after=Inf)) %>% 
                   mutate(released_test = fct_recode(released_test,
                                                             "No post-\nflight test"  = "No post-flight test",
                                                             "One post-\nflight test" = "One post-flight test")) %>% 
                   #filter(pre_board_screening!="No pre-flight test")%>% 
                   filter(stringency!="high") ,
                 aes(x=time_in_iso,
                     ymin=`0.025`,
                     ymax=`0.975`,
                     colour=released_test),
                 position = position_dodge2(width=0.9),
                 size=4,
                 alpha=0.7)+
  geom_linerange(data=plot_data %>% mutate(stringency=fct_recode(stringency,"max"="maximum")) %>% 
                   mutate(stringency=fct_relevel(as.factor(stringency),"low","moderate","high","max",after=Inf)) %>% 
                   mutate(released_test = fct_recode(released_test,
                                                             "No post-\nflight test"  = "No post-flight test",
                                                             "One post-\nflight test" = "One post-flight test")) %>% 
                   #filter(pre_board_screening!="No pre-flight test") %>% 
                   filter(stringency!="high") ,
                 aes(x=time_in_iso,
                     ymin=`0.25`,
                     ymax=`0.75`,
                     colour=released_test),
                 position = position_dodge2(width=0.9),
                 size=4)+
  scale_colour_manual(name="",
                      values=lshtm_greens)+
  ggnewscale::new_scale_colour()+
  geom_linerange(data=plot_data %>% 
                   filter(stringency=="high") %>% mutate(stringency=fct_recode(stringency,"max"="maximum")) %>% 
                   mutate(stringency=fct_relevel(as.factor(stringency),"low","moderate","high","max",after=Inf)) %>% 
                   mutate(released_test = fct_recode(released_test,
                                                             "No post-\nflight test"  = "No post-flight test",
                                                             "One post-\nflight test" = "One post-flight test")) %>% 
                  # filter(pre_board_screening!="No pre-flight test") %>% 
                   left_join(bivariate_color_scale),
                 aes(x=time_in_iso,
                     ymin=`0.025`,
                     ymax=`0.975`,
                     colour=colour),
                 position = position_dodge2(width=0.9),
                 size=4,
                 alpha=0.7)+
  geom_linerange(data=plot_data %>% 
                   filter(stringency=="high") %>% mutate(stringency=fct_recode(stringency,"max"="maximum")) %>% 
                   mutate(stringency=fct_relevel(as.factor(stringency),"low","moderate","high","max",after=Inf)) %>% 
                   mutate(released_test = fct_recode(released_test,
                                                             "No post-\nflight test"  = "No post-flight test",
                                                             "One post-\nflight test" = "One post-flight test")) %>% 
                   #filter(pre_board_screening!="No pre-flight test") %>% 
                   left_join(bivariate_color_scale),
                 aes(x=time_in_iso,
                     ymin=`0.25`,
                     ymax=`0.75`,
                     colour=colour),
                 position = position_dodge2(width=0.9),
                 size=4)+
  scale_colour_identity(guide=FALSE)+
  scale_x_continuous(breaks=breaks_width(1))+
  scale_y_continuous(expand = expansion(mult = c(0,.1), add = c(0,0) )
  )+
  facet_nested(country~pre_board_screening+stringency,scales="free_x",
               nest_line = TRUE,
               switch = NULL,labeller = labeller(stringency=capitalize)) +
  force_panelsizes(cols=c(2,6,8,2))+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA),
        text = element_text(size=20),
        strip.placement = "outside",
        panel.grid.minor.x = element_blank(),
        legend.position = c(0.82,0.95))+
  xlab("Days in quarantine")+
  ylab("Number of infectious persons released per week")

fig3 <- ggdraw()+
  draw_plot(pre_board_plot,0,0,1,1)+
  draw_plot(legend,0.78,0.78,0.25,0.15)

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig3.",.x),
              plot=fig3,
              width = 600, height = 300, units="mm",
              dpi = 280))


