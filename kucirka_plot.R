
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
                 aes(y = ..density..,
                     fill=type),
                 alpha=0.8,
                 colour = "white") +
  scale_fill_manual(values=lshtm_greens,guide=F)+
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


fig1 <- kucirka_plot / distributions_plot + plot_layout(ncol = 1, guide = "collect",heights = c(1,2)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(legend.position = "bottom"))

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/fig1_.",.x),
              plot     = fig1,
              width    = 210,
              height   = 210,
              units    = "mm",
              dpi      = 300))
### end distributions
