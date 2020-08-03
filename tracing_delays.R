### Test and Trace delays
source("packages.R")

#Delay from positive test to getting info on close contacts from index
getting_contact_info <- 
  readxl::read_excel("data/NHS_T_T_timeseries_Template_Output_Week_8.xlsx",
                     sheet = "Table_9", 
                     skip = 2) %>% 
  dplyr::rename(desc = `...1`) %>%
  dplyr::filter(grepl(pattern = "Number", x = desc)) %>%
  select(desc, contains("Week")) %>%
  select(1, n = ncol(.)) %>%
  mutate(t=c(0.5,1.5,2.5,3.5)) %>% 
  mutate(n=as.numeric(n)) %>% 
  uncount(n)

#Delay from getting info to tracing contacts
tracing_delay <-
  readxl::read_excel("data/NHS_T_T_timeseries_Template_Output_Week_8.xlsx",
                     sheet = "Table_12",
                     skip = 2) %>% 
  dplyr::rename(desc = `...1`) %>%
  dplyr::filter(grepl(pattern = "Number", x = desc)) %>%
  select(desc, contains("Week")) %>%
  select(1, n = ncol(.)) %>%
  mutate(t=c(0.5,1.5,2.5,3.5)) %>% 
  mutate(n=as.numeric(n)) %>% 
  uncount(n)



