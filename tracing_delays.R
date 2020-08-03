### Test and Trace delays
source("packages.R")

#Delay from positive test to getting info on close contacts from index
getting_contact_info <- readxl::read_excel("data/NHS_T_T_timeseries_Template_Output_Week_8.xlsx",sheet = "Table_9") %>% 
  slice(-1:-3) %>% 
  slice(1:4) %>% 
  select(c(1,18)) %>% 
  rename("desc"=1,
         "n"=2) %>% 
  mutate(t=c(0.5,1.5,2.5,3.5)) %>% 
  mutate(n=as.numeric(n)) %>% 
  uncount(n)

#Delay from getting info to tracing contacts
tracing_delay <- readxl::read_excel("data/NHS_T_T_timeseries_Template_Output_Week_8.xlsx",sheet = "Table_12") %>% 
  slice(-1:-3) %>% 
  slice(1:4) %>% 
  select(c(1,18)) %>% 
  rename("desc"=1,
         "n"=2) %>% 
  mutate(t=c(0.5,1.5,2.5,3.5)) %>% 
  mutate(n=as.numeric(n)) %>% 
  uncount(n)



