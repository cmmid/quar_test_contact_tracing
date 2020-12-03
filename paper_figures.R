## curves for paper

set.seed(145)
trajectories <- make_trajectories(n_cases = 10)

trajectories_to_plot <- trajectories$models %>%
  filter(type == "symptomatic") %>%
  mutate(pred = map(.x = m, 
                    ~data.frame(x = seq(0, 25, length.out = 1001)) %>%
                      mutate(y = .x(x)))) %>%
  unnest(pred) %>%
  mutate(val = base::cut(y, 
                         breaks = c(-Inf, 27, 30, 35, Inf),
                         labels = c("0.95", "0.75", "0.3", "0")),
         LFA = parse_number(as.character(val))) %>%
  mutate(PCR = as.integer(y < 35)) %>%
  gather(key, value, PCR, LFA)


mask <- data.frame(min = c(-Inf, 27, 30, 35),
                   max = c(27, 30, 35, Inf),
                   val = c(.95, .75, .3, 0),
                   key = "LFA") %>%
  bind_rows(data.frame(min = c(-Inf, 35),
                       max = c(35, Inf),
                       val = c(1, 0),
                       key = "PCR"))

trajectories_plot <- trajectories_to_plot %>%
  ggplot(data = ., aes(x = x, y = y)) +
  annotate(geom = "rect",
           ymax = 30, xmin = 0, xmax = 25, ymin = -Inf,
           color = NA,
           fill = "black", alpha = 0.1) +
  geom_line(aes(group = idx,
                color = factor(value)),
            alpha = 0.75) +
  theme_bw() +
  facet_grid(. ~ key) +
  geom_hline(data = mask,
             aes(yintercept = max),
             lty = 2) +
  scale_color_manual(values = c(muted("red", l = 80), 
                                RColorBrewer::brewer.pal(4, "Purples")[-1],
                                muted("blue", l = 50, c = 150)),
                     name   = "Probability of detection") +
  ylab(expression(C[t]~value)) +
  xlab("Time since exposure (days)") +
  plotting_theme +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2) )) +
  geom_text(x = 25, y = 29,
            hjust = 1,  size = 3,
            label = "Infectious") +
  geom_text(x = 25, y = 31,
            hjust = 1, size = 3,
            label = "Non-infectious")

save_plot(plot = trajectories_plot, 
          device="png",prefix = "trajectories", width = 210,height=105)


curves <- read_rds("data/matched_curves.rds")


curves %>% filter(rowid < 10) %>%
  ggplot(data= ., aes(x = days_since_infection, y = value)) +
  geom_line(aes(color = factor(assay))) + facet_wrap(~rowid) +
  plotting_theme +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  xlab("Time since exposure (days)") +
  ylab("Probability of detection") +
  ylim(c(0,1)) +
  scale_color_
