X <- ls(pattern = "^result_") %>%
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

delay_pred <- data.frame(x = seq(0, 4, by = 0.01))

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
  mutate(index_result_delay = 
           time_to_event(n = n(),
                         mean = P_r[["mean"]],
                         var  = P_r[["var"]])) %>% 
  #sample contact info delay
  mutate(contact_info_delay = 
           time_to_event(n = n(),
                         mean = P_c[["mean"]],
                         var  = P_c[["var"]])) %>% 
  #sample tracing delay
  mutate(tracing_delay      =
           time_to_event(n = n(),
                         mean = P_t[["mean"]],
                         var  = P_t[["var"]])) %>%
  mutate(total = index_result_delay + contact_info_delay + tracing_delay)

total_density <- 
  pull(delay_total, total) %>%
  {density(log(.))} %>%
  {data.frame(x = exp(.$x),
              y = .$y, variable = "Total")}

total_quantiles <- delay_total %>%
  select(total) %>%
  mutate(variable = "Total") %>%
  nest(data = total) %>%
  mutate(Q = map(data, ~quantile(.x$total, probs) %>% set_names(., probs))) %>%
  unnest_wider(Q) %>%
  mutate(mean = map_dbl(data, ~mean(.x$total)),
         var  = map_dbl(data, ~var(.x$total)))

delay_distributions <- 
  list(
    `Testing` = index_result_delay,
    `Sourcing`= getting_contact_info,
    `Tracing` = tracing_delay) %>%
  map_df(~data.frame(delay_to_gamma(.x)) %>%
           mutate(n = nrow(.x)) %>%
           bind_cols(mv2gamma(mean = .$mean, var = .$var)),
         .id = "variable")  %>%
  mutate(variable = fct_inorder(factor(variable))) %>%
  nest(data = c(scale, shape, rate)) %>%
  mutate(Delay = map(data, ~mutate(delay_pred,
                                   y = dgamma(x = x, 
                                              shape = .x$shape,
                                              rate  = .x$rate)))) %>%
  unnest(Delay) %>%
  bind_rows(total_density) %>%
  mutate(variable = fct_inorder(factor(variable)))

delay_quantiles <- 
  delay_distributions %>% 
  filter(x == 0,
         variable != "Total") %>%
  mutate(Q = map(data, 
                 ~qgamma(p = probs, 
                         shape = .x$shape,
                         rate = .x$rate) %>%
                   set_names(., probs)) 
  ) %>%
  unnest_wider(Q) %>%
  bind_rows(total_quantiles) %>%
  mutate(variable = fct_inorder(factor(variable)))


delay_histograms <- 
  list(
    `Testing` = index_result_delay,
    `Sourcing`= getting_contact_info,
    `Tracing` = tracing_delay) %>%
  bind_rows(.id = "variable")  %>%
  mutate(variable = fct_inorder(factor(variable))) %>%
  count(variable, t) %>%
  group_by(variable) %>%
  mutate(y = n/sum(n)) %>%
  rename(x = t)

delay_labeller <- function(x){
  fct_recode(x,
             "Return of index case's test" = "Testing",
             "Sourcing of contacts"        = "Sourcing",
             "Tracing of contacts"         = "Tracing",
             "Total"                       = "Total")
}

delay_plot <- ggplot(data = delay_distributions, aes(x=x, y=y)) + 
  geom_col(data = delay_histograms,
           fill = lshtm_greens[2],
           color = NA, width = 1) +
  geom_line(color = lshtm_greens[1]) +
  geom_segment(data = delay_quantiles,
               aes(x = `0.025`, xend = `0.975`,
                   y = 0, yend = 0),
               color = lshtm_greens[1],
               size = 1) +
  geom_point(data = delay_quantiles,
             aes(x = `0.5`, y = 0),
             pch = 21,
             fill = "white",
             color = lshtm_greens[1],
             size = 1) +
  facet_wrap( ~variable,
              nrow = 1,
              labeller = labeller(variable = delay_labeller)) +
  theme_bw() +
  xlab("Time to completion (days)") +
  ylab("Density") +
  coord_cartesian(ylim = c(0, 1.5),
                  xlim = c(0, 6))

save_plot(delay_plot, 
          prefix = "delays",
          base = "plot", 
          width = 210, 
          height = 40)




