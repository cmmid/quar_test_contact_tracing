X <- ls(pattern = "^P") %>%
  set_names(.,.) %>%
  map(get) %>%
  map_df(~data.frame(delay_to_gamma(.x)) %>%
           mutate(n = nrow(.x)) %>%
           bind_cols(mv2gamma(mean = .$mean, var = .$var)),
         .id = "source")

X_for_plot <- ls(pattern = "^result_") %>%
  set_names(.,.) %>%
  map_df(get, .id = "source") %>%
  count(source, t) %>%
  group_by(source) %>%
  mutate(y = n/sum(n)) %>%
  rename(x = t)

# what sort of delays do we see when breaking down test results by source?
Table_1a_summaries <- 
  X %>%
  nest(data = c(rate, shape)) %>%
  mutate(Q = map(.x = data, ~with(.x, qgamma(p = probs,
                                             shape = shape,
                                             rate  = rate)) %>%
                   set_names(probs))) %>%
  unnest_wider(Q) %>%
  unnest_wider(data)

delay_pred <- data.frame(x = seq(0, 10, by = 0.01))

source_labeller <- function(x){
  sub(pattern = "result_delay_", replacement = "", x = x) %>%
    capitalize
}

Table_1a_summaries %>%
  select(source, mean, var, scale, rate, shape) %>%
  nest(data = c(scale, shape, rate)) %>%
  mutate(Delay = map(data, ~mutate(delay_pred,
                                   y = dgamma(x = x, 
                                              shape = .x$shape,
                                              rate  = .x$rate)))) %>%
  unnest(Delay) %>%
  ggplot(data = ., aes(x=x, y=y)) + 
  geom_col(data = X_for_plot,
           fill = lshtm_greens[2], color = NA, width = 1) +
  geom_line(color = lshtm_greens[1]) +
  facet_wrap( ~source, labeller = labeller(source = source_labeller)) +
  theme_bw() +
  xlab("Delay to testing result") +
  ylab("Density")



## distributions for paper

## total



delay_total <- data.frame(sim = 1:10000) %>%
  mutate(Testing = 
           time_to_event(n = n(),
                         mean = P_r[["mean"]],
                         var  = P_r[["var"]])) %>% 
  #sample contact info delay
  mutate(Sourcing = 
           time_to_event(n = n(),
                         mean = P_c[["mean"]],
                         var  = P_c[["var"]])) %>% 
  #sample tracing delay
  mutate(Tracing      =
           time_to_event(n = n(),
                         mean = P_t[["mean"]],
                         var  = P_t[["var"]])) %>%
  mutate(Total = Testing + Sourcing + Tracing)

total_density <- 
  pull(delay_total, Total) %>%
  {density(log(.))} %>%
  {data.frame(x = exp(.$x),
              y = .$y, Event = "Total")}


delay_histograms <- 
  list(
    `Testing` = index_result_delay,
    `Sourcing`= getting_contact_info,
    `Tracing` = tracing_delay) %>%
  bind_rows(.id = "Event")  %>%
  mutate(Event = fct_inorder(factor(Event))) %>%
  count(Event, t) %>%
  group_by(Event) %>%
  mutate(y = n/sum(n)) %>%
  rename(x = t)

delay_labeller <- function(x){
  fct_recode(x,
             "Return of index case's test" = "Testing",
             "Sourcing of contacts"        = "Sourcing",
             "Tracing of contacts"         = "Tracing",
             "Total"                       = "Total")
}

delay_distributions <-
  list(
    Testing  = P_r,
    Sourcing = P_c,
    Tracing  = P_t) %>%
  map(~mv2gamma(.x$mean, .x$var)) %>%
  map_df(.x = .,
         ~mutate(delay_pred,
                 y = dgamma(x = x,
                            shape = .x$shape,
                            rate  = .x$rate)),
         .id = "Event") %>%
  bind_rows(total_density) %>%
  mutate(Event = fct_inorder(factor(Event)))

delay_plot <- ggplot(data = delay_distributions, aes(x=x, y=y)) + 
  geom_col(data = delay_histograms,
           fill = lshtm_greens[2],
           color = NA, width = 1) +
  geom_line(color = lshtm_greens[1]) +
  geom_segment(data = delay_summary,
               aes(x = `2.5%`, xend = `97.5%`,
                   y = 0, yend = 0),
               color = lshtm_greens[1],
               size = 1) +
  geom_point(data = delay_summary,
             aes(x = `50%`, y = 0),
             pch = 21,
             fill = "white",
             color = lshtm_greens[1],
             size = 1) +
  facet_wrap( ~Event,
              nrow = 2,
              labeller = labeller(Event = delay_labeller)) +
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA)) +
  xlab("Time to completion (days)") +
  ylab("Density") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits = c(0,12),
                     breaks = pretty_breaks(4))


save_plot(delay_plot, dpi = 320, device = "png",
          prefix = "delays",
          base = "plot", 
          width = 210, 
          height = 140)

ashcroft_plot <- time_to_event_lnorm(n = 1e4, 
                                     inc_parms$mu_inc,
                                     inc_parms$sigma_inc) %>%
  {data.frame(index_onset_t = .)} %>%
  mutate(sec_exposed_t = index_onset_t - infect_shift + 
           rtgamma(n     = n(),
                   a     = infect_shift,
                   b     = Inf, 
                   shape = infect_shape,
                   rate  = infect_rate) 
  ) %>%
  ggplot(data = ., aes(x = sec_exposed_t)) +
  geom_histogram(binwidth = 1, center = 0.5,
                 aes(y = ..density..),
                 fill = lshtm_greens[2],
                 color = lshtm_greens[1],
                 alpha = 0.7)  +
  theme_minimal() + 
  theme(panel.border=element_rect(fill=NA)) +
  xlab("Time from exposure (days)") +
  ylab("Density") +
  #coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(limits = c(0,20),
                     breaks = pretty_breaks(4))

save_plot(ashcroft_plot, dpi = 320, device = "png",
          prefix = "ashcroft",
          base = "plot", 
          width = 210, 
          height = 140)
