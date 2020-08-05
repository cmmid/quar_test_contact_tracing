source("packages.R")
source("utils.R")
source("tracing_delays.R")
source("wolfel.R")
source("he.R")

input <- 
  tibble(pathogen = "SARS-CoV-2") %>%
  bind_cols(., list(
    `low` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = 0,
               second_test_delay = NA), 
    `moderate` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = c(3,5,7),
               second_test_delay = NA),
    `high` = 
      crossing(screening         = TRUE,
               first_test_delay  = c(0:3),
               second_test_delay = c(2,4,6)),
    `maximum` = 
      crossing(screening         = c(TRUE,FALSE),
               first_test_delay  = 14,
               second_test_delay = NA)) %>%
      bind_rows(.id = "stringency")) %>% 
  mutate(scenario=row_number()) %>%
  mutate(max_mqp             = 14,
         post_symptom_window =  7,
         results_delay       =  1
  )
baseline_low <- data.frame(
  screening = FALSE,
  first_test_delay      = 0,
  second_test_delay     = NA,
  stringency            = "low"
)

baseline_max <- data.frame(
  screening = FALSE,
  first_test_delay      = 14,
  second_test_delay     = NA,
  stringency            = "maximum"
)

results_df <- tibble(index_test_delay=c(0:4)) %>% 
     mutate(results =map(.f=run_analysis,.x=index_test_delay,
                         n_sims          = 100,
                         n_sec_cases     = 1000, # this shouldn't matter. just needs to be Big Enough
                         n_ind           = 10000,
                         seed            = 145,
                        index_result_delay = index_result_delay,
                        contact_info_delay = getting_contact_info,
                        tracing_delay      = tracing_delay,
                        asymp_parms        = asymp_fraction))

# need to spit this out to file rather than default graphics
results_df %<>% 
  mutate(rr_plot_dat = pmap(.l=list(released_times=results),
                                  .f=run_rr_analysis,
                            
                         main_scenarios = main_scenarios,
                         baseline_scenario = baseline_max,
                         text_size = 2.5,
                         faceting = ~ stringency))
 
ylabA <- sprintf("Ratio of infectious persons released in comparison to\n%s stringency, %i day quarantine, %s scenario",
                 baseline_scenario$stringency,
                 with(baseline_scenario,
                      first_test_delay + screening + 
                        ifelse(is.na(second_test_delay), 0, second_test_delay)),"no testing")

results_fig <- results_df %>% 
  select(-results) %>%  
  unnest(rr_plot_dat) %>% 
  plot_data(input = input) %>% 
  make_release_figure(text_size = 2.5,
                      xlab      = "x",
                      ylab      = "Ratio of infectious persons released in comparison to\n%s stringency, %i day quarantine, %s scenario",, 
                      hline="dashed",
                      log_scale = TRUE,
                      faceting  = index_test_delay~stringency)

list("png", "pdf") %>%
  map(~ggsave(filename = paste0("results/rr_figs_baseline_",
                                "max",".",.x),
              plot=results_fig,
              width = 297, 
              height = 297*1.5,#120*nrow(distinct(ungroup(n_risk_ratios),
                        #                 !!lhs(faceting))), 
              units="mm",
              dpi = 400,
              device = ifelse(.x=="pdf",cairo_pdf,
                              "png")))



