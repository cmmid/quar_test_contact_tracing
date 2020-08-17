faceting <- index_test_delay + delay_scaling ~ stringency

results_df <- results %>% 
  map(~mutate(.x, infectivity_total   = (infectivity_post + infectivity_pre),
         infectivity_averted = 1 - infectivity_total))

results_infectivity <- 
  results_df %>%
  make_days_plots(.,
      faceting = faceting,
      y_labels =
        c("infectivity_pre" = 
            "Transmission potential prior to secondary cases being traced",
          "infectivity_post" =
            "Transmission potential after release of secondary cases",
          "infectivity_averted" = 
            "Transmission potential averted as a result of quarantine and testing of secondary cases"
        ),
      base = "all",
      sum = F)

c("infectivity_post" =
    "Transmission potential after release of secondary cases",
  "infectivity_averted" = 
    "Transmission potential averted as a result of quarantine and testing of secondary cases"
) %>%
  map2(.x = ., .y = names(.),
       .f = ~set_names(.x, .y)) %>%
  map(
    ~make_days_plots(results_df,
                     input, 
                     faceting = faceting,
                     y_labels = .x,
                     base = names(.x),
                     sum = F))




results_df %>%
  make_days_plots(input, 
                  faceting = faceting,
                  y_labels = c("infectivity_total" = 
                                 "Transmission potential in community compared to no quarantine or testing of secondary cases"),
                  base = "total",
                  sum = F)


results_df %>%
  make_days_plots(input, 
                  faceting = faceting,
                  y_labels = c("infectivity_averted" = 
                                 "Transmission potential averted as a result of quarantine and testing of secondary cases"),
                  base = "averted",
                  sum = F)



results_infectivity_df <-
  results_infectivity %>% 
  map(show_results, reduction = FALSE) %>%
  map_df(bind_rows, .id = "Measure")

# transmission potential ahead of tracing
results_infectivity_df %>%
  filter(stringency == "low", 
         screening == FALSE,
         Measure == "pre") %>%
  select(one_of(all.vars(faceting)),
         contains("%")) %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~percent(round(., 2)))

# how much quarantine time can we shave off if we improve testing delays?
c(2,3) %>%
  set_names(., paste(., "days until index case's test")) %>%
  map(~filter(results_infectivity_df, Measure == "averted") %>%
        filter((stringency == "maximum" & screening == FALSE & index_test_delay == .x) |
                 (stringency == "moderate" & index_test_delay == .x - 1)) %>%
        select(-Measure, -delays, -screening) %>%
        mutate_at(.vars = vars(contains("%")),
                  .funs = ~percent(round(., 2))))

# if we can't reduce delays, can we do double testing?
c(1, 2, 3) %>%
  set_names(., paste(., "days until index case's test")) %>%
  map(~filter(results_infectivity_df, Measure == "averted" & index_test_delay == .x) %>%
        filter((stringency == "maximum" & screening == FALSE) |
                 (stringency == "moderate" )) %>%
        select(-Measure, -screening) %>%
        mutate_at(.vars = vars(contains("%")),
                  .funs = ~percent(round(., 2))))

# baseline at Maximum
results_infectivity_df %>%
  filter(index_test_delay == 2,
         stringency == "maximum",
         Measure == "averted") 

results_infectivity_df %>%
  filter(screening == FALSE,
         stringency == "maximum",
         Measure == "averted")  %>%
  select(index_test_delay, time_in_iso, contains("%"))

# Table 2: max and 10 days (with one or two tests)
results_infectivity_df %>%
  filter(screening == TRUE,
         time_in_iso == 10,
         Measure == "averted")  %>%
  select(stringency, index_test_delay,delays, time_in_iso, contains("%"))  %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(.,2)) %>%
  select(-time_in_iso)

# baseline with test on day 0
results_infectivity_df %>%
  filter(index_test_delay == 2,
         stringency == "low",
         Measure == "averted") %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(.,2)) %>%
  select(-time_in_iso)

# effect of a test in low
results_infectivity$averted %>%
  filter(stringency == "moderate") %>%
  show_results(reduction = TRUE)

# effect of asymptomatics

results_infectivity_type <- 
  results_df %>%
  make_days_plots(input, 
                  faceting = index_test_delay  ~ stringency + type,
                  y_labels =
                    c("infectivity_averted" = 
                        "Transmission potential averted as a result of quarantine and testing of secondary cases"),
                  base = "averted_type_",
                  sum = F)

# Table 3... Table 2 + type
results_infectivity_type$averted %>%
  split(.$type) %>%
  map(~filter(.x,time_in_iso == 14 | time_in_iso == 1 | 
                (time_in_iso == 10 & stringency == "moderate")) %>%
        show_results(reduction = F)) %>%
  map_df(bind_rows, .id = "type") %>%
  mutate_at(.vars = vars(contains("%")),
            .funs = ~round(100 * .)) %>%
  arrange(desc(time_in_iso), index_test_delay, desc(type)) %>%
  mutate(`PCR tests` = 
           ifelse(!screening, 
                  "no test",
                  paste("day", sub(pattern = " & NA", 
                                   replacement = "", 
                                   x = delays)))) %>%
  mutate(`PCR tests` = sprintf("%s (%s)", time_in_iso, `PCR tests`)) %>%
  select(-screening, -delays, -time_in_iso) %>%
  group_by_at(.vars = vars(-contains("%"))) %>%
  transmute(Median = sprintf("%.f%%", `50%`),
            `50% CI (IQR)`  = sprintf("(%.f%%, %.f%%)", `25%`, `75%`),
            `95% CI`  = sprintf("(%.f%%, %.f%%)", `2.5%`, `97.5%`)) %>%
  write_csv(., "results/Table_3.csv")
