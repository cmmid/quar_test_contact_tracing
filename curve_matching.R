## matching PCR and LFA curves

PCR_curves <- read_csv("data/posterior_samples_ct_threshold_37.csv")
LFA_curves <- read_csv("data/posterior_samples_ct_threshold_28.csv")

list(LFA = LFA_curves,
     PCR = PCR_curves) %>%
  map_df(~filter(.x, iter <= 10) %>%
           select(-X1), .id = "Assay") %>%
  pivot_wider(names_from = "Assay", values_from = "value") %>%
  ggplot(data = ., aes(x=PCR, y = LFA)) +
  geom_path(aes(group = iter)) +
  facet_wrap(~iter)


standardise <- function(x){
  (x - mean(x))/sd(x)
}

my_dist <- function(x, y){
  #  browser()
  X <- unlist(x[["diff_r"]]  - y[["diff_r"]])
  Y <- unlist(x[["value_r"]] - y[["value_r"]])
  
  D <- sqrt(X^2 + Y^2)
  D_min <- D == min(D)
  
  # to break ties, sample at random
  sample_n(select(filter(y, D_min), iter_pcr = iter), size = 1)
}

list(LFA = LFA_curves,
     PCR = PCR_curves) %>%
  map(~group_by(.x, iter) %>%
        filter(value == max(value)) %>%
        ungroup) %>%
  map(~head(.x, 100)) %>%
  map(~mutate(.x, 
              value_s = standardise(value),
              diff_s  = standardise(diff),
              value_r = rank(value),
              diff_r  = rank(diff))) %>%
  {bind_cols(.[[1]],
             bind_rows(lapply(X = group_split(rowwise(.[[1]])),
                              FUN = function(x){
                                my_dist(x, .[[2]])
                              })
             )) 
  } %>%
  rename(iter_lfa = iter) -> mapping

mapping_plot <- 
  mapping %>%
  select(iter_lfa, iter_pcr) %>%
  tibble::rowid_to_column(.) %>%
  nest(data = -c(rowid, iter_lfa)) %>%
  inner_join(LFA_curves, by = c("iter_lfa" = "iter")) %>%
  unnest(data) %>%
  select(-X1) %>%
  rename(value_lfa = value) %>%
  left_join(PCR_curves, by = c("iter_pcr" = "iter", "diff")) %>%
  rename(value_pcr = value) %>%
  gather(key, value, value_lfa, value_pcr) %>%
  ggplot(data = ., aes(x = diff, y = value)) +
  geom_line(aes(color = key)) +
  facet_wrap(~rowid) +
  theme_void() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

ggsave(filename = "A.pdf", plot = mapping_plot)


# now that we have iter_lfa and iter_pcr we have a map to get the nearest PCR for each LFA curve

