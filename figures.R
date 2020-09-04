# main figures

faceting      <- index_test_delay + delay_scaling + waning + stringency  ~ .
faceting_wide <- index_test_delay + delay_scaling ~ waning + stringency

infectivity_labels <-
  c("infectivity_post" =
      "Transmission potential of secondary cases \nafter release",
    "infectivity_averted" = 
      "Transmission potential of secondary cases \naverted as a result of quarantine and testing",
    # "infectivity_quar" = 
    #   "Transmission potential in community due to imperfect quarantine adherence",
    "infectivity_pre" =
      "Transmission potential of secondary cases \nprior to being traced",
    "infectivity_total" = 
      "Transmission potential of secondary cases \nin community compared to no quarantine or testing"
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
  {.[filter(input,index_test_delay == 2,
           waning == "waning_constant") %>% pull(scenario)]} %>%
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
  .[2] %>%
  map2(.x = ., .y = names(.),
       .f = ~set_names(.x, .y)) %>%
  map(
    ~make_days_plots(x = get(results_name),
                      input, 
                     plot = T,
                     faceting = faceting_wide,
                     y_labels = .x,
                     dir = results_name,
                     base = paste(results_name, names(.x),sep="_"),
                     sum = F))


