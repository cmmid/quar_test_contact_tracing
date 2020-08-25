# main figures

faceting <- index_test_delay + delay_scaling + waning ~ stringency

infectivity_labels <-
  c("infectivity_post" =
      "Transmission potential after release of secondary cases",
    "infectivity_averted" = 
      "Transmission potential averted as a result of quarantine and testing of secondary cases",
    "infectivity_quar" = 
      "Transmission potential in community due to imperfect quarantine adherence",
    "infectivity_pre" =
      "Transmission potential prior to secondary cases being traced",
    "infectivity_total" = 
      "Transmission potential in community compared to no quarantine or testing of secondary cases"
  )

#results <- readRDS("results/results.RDS") 

# results %<>%
#   map(~mutate(.x,
#               infectivity_total = (infectivity_post + 
#                                          infectivity_pre),
#               infectivity_averted = 1 - infectivity_total))

#results_df <- results

results_infectivity <- 
  get(results_name) %>%
  make_days_plots(.,
      faceting = faceting,
      y_labels = grep(value = T, pattern = "prior",
                      x = infectivity_labels, invert = T),
      dir = results_name,
      base = "all",# all
      sum = F)




# results_waning %>%
#   make_days_plots(.,
#                   faceting = faceting,
#                   y_labels = infectivity_labels["infectivity_quar"],
#                   base = "waning_quar_cons", # all
#                   sum = F)

infectivity_labels %>%
  map2(.x = ., .y = names(.),
       .f = ~set_names(.x, .y)) %>%
  map(
    ~make_days_plots(get(results_name),
                     input, 
                     faceting = faceting,
                     y_labels = .x,
                     dir = results_name,
                     base = paste(results_name, names(.x),sep="_"),
                     sum = F))


