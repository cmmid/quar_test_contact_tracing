## curves for paper

set.seed(145)
trajectories <- make_trajectories(n_cases = 10)

trajectories_to_plot <- trajectories$models %>%
  filter(type == "symptomatic") %>%
  mutate(pred = map(.x = m, 
                    ~data.frame(x = seq(0, 25, length.out = 1001)) %>%
                      mutate(y = .x(x)))) %>%
  unnest(pred) %>%
  mutate(LFA = base::cut(y, 
                         breaks = c(-Inf,20,25,30,35,Inf),
                         labels=c("82.4%","54.5%","8.3%","5.3%","0%"))) %>%
  mutate(PCR = percent(as.integer(y < 35))) %>%
  gather(key, value, PCR, LFA) %>% 
  mutate(value=fct_rev(fct_relevel(value,"100%","82.4%","54.5%","8.3%","5.3%","0%")))

mask <- data.frame(min = c(-Inf, 20, 25, 30, 35),
                   max = c(20, 25, 30, 35, Inf),
                   val = c(.824, .545, .083,0.053,0),
                   key = "LFA") %>%
  bind_rows(data.frame(min = c(-Inf, 35),
                       max = c(35, Inf),
                       val = c(1, 0),
                       key = "PCR"))

trajectories_to_plot %>%
  ggplot(data = ., aes(x = x, y = y)) +
  annotate(geom = "rect",
           ymax = Inf, xmin = -Inf, xmax = Inf, ymin = 30,
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
  scale_color_manual(values = viridis::magma(n=10,direction=-1)[c(2:8)],
                     name   = "Probability of detection") +
  ylab("Ct value") +
  xlab("Time since exposure (days)") +
  plotting_theme +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2),nrow = 1)) +
  geom_text(x = 25, y = 29,
            hjust = 1,
            label = "Infectious") +
  geom_text(x = 25, y = 31,
            hjust = 1, 
            label = "Non-infectious")

save_plot(dpi=400,
          device="png",prefix = "trajectories", width = 210,height=105)
