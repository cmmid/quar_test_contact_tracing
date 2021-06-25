# Load required packages and utility scripts
source("scripts/packages.R")
source(here("scripts","utils.R"))
source(here("scripts","plot_functions.R"))
source(here("scripts","parameters.R"))


most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_act_a.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

results_df <- read.fst(most_recent_file)

#results table
results_df %>% 
  #filter(!is.infinite(inf_start) & !is.infinite(inf_end)) %>%
  mutate(inf_period = 
           ifelse(is.infinite(inf_start) & is.infinite(inf_end),0, inf_end-inf_start),
         test_num = 
           case_when(
             stringency ==        1                     ~ 0,
             stringency ==        2                     ~ 1,
             stringency ==        3 & test_no != "None" ~ parse_number(test_no) + 1,
             stringency ==        3 & test_no == "None" ~ 1,
             stringency %in% c(4,5) & test_no != "None" ~ parse_number(test_no),
             stringency %in% c(4,5) & test_no == "None" ~ n_tests,
             TRUE ~ NA_real_) # just in case!
  ) %>% 
  group_by(ind_idx, stringency, assay,
           test_exit_self_iso, test_to_tracing,  
           quar_dur, n_tests, 
           iso_dur) %>% 
  summarise(all        = sum(inf_period),
            in_iso     = sum(symp_overlap+test_overlap),
            in_quar    = sum(quar_overlap),
            in_comm    = all-(in_iso+in_quar),
            tests_used = sum(test_num)
  ) %>% 
  pivot_longer(cols=c(all,in_quar,in_iso,in_comm,tests_used)) %>% 
  group_by_at(.vars = vars(-value, -ind_idx)) %>%
  
  #group_by(stringency,assay,quar_dur,test_exit_self_iso,n_tests,test_to_tracing,iso_dur,name) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$value,
                              probs = c(0.025,0.5,0.975)))) %>%
  unnest_wider(Q) %>% 
  mutate(stringency = factor(stringency)) %>% 
  crossing(contacts = 10000,
           prevalence = c(0.01, 0.1, 0.5)) %>% 
  mutate(n_inf = contacts*prevalence,
         across(c(`2.5%`:`97.5%`),function(x){(round(x*n_inf/10,0))}),
         uninf_test_num=(contacts-n_inf)*case_when(stringency ==        1 ~ 0,
                                                   stringency ==        2 ~ 1,
                                                   stringency ==        3 ~ 1,
                                                   stringency %in% c(4,5) ~ n_tests)
  ) %>% 
  select(-data) %>% 
  pivot_wider(values_from = c(`2.5%`:`97.5%`),
              names_from  = name,
              names_glue  = "{name}_{.value}",
              names_sort  = F) %>% 
  mutate(across(c(`tests_used_50%`,
                  `tests_used_2.5%`,
                  `tests_used_97.5%`),
                function(x){x+uninf_test_num})) %>%
  arrange(stringency, prevalence, test_to_tracing, quar_dur, n_tests, iso_dur) %>% 
  ungroup() %>% 
  relocate(
    `all_50%`,`all_2.5%`,`all_97.5%`,
    `in_quar_50%`,`in_quar_2.5%`,`in_quar_97.5%`,
    `in_iso_50%`,`in_iso_2.5%`,`in_iso_97.5%`,
    `in_comm_50%`,`in_comm_2.5%`,`in_comm_97.5%`,
    `tests_used_50%`,`tests_used_2.5%`,`tests_used_97.5%`,
    everything()) %>% 
  write.csv("results/release_results.csv")
