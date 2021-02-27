# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")

results_name <- "results_df"

most_recent_file <- file.info(list.files("results/", full.names = T)) %>% 
  as.data.frame() %>% 
  rownames_to_column()%>% 
  filter(str_detect(rowname,"_main.fst")) %>% 
  slice_max(mtime) %>% 
  pull(rowname)

assign(results_name,read.fst(most_recent_file))

select(-c(adherence_quar,adherence_iso)) %>% 
  crossing(adherence_iso=seq(0,1,by=0.1),adherence_quar=seq(0,1,by=0.1)) %>% 
  mutate(adhering_quar=rbinom(n=n(),size = 1,prob = adherence_quar),
         adhering_iso=rbinom(n=n(),size = 1,prob = adherence_iso)) %>% 
  group_by(adherence_quar,adherence_iso,test_to_tracing) %>%
  group_split() 

plan(multisession, workers = 7)
adherence_results <- adherence_list %>% 
  future_map(.x=.,.f=calc_overlap,.progress = T)

a <- get(results_name) %>% 
  mutate(adherence_quar=as.factor(adherence_quar)) %>% 
  mutate(lft_delivery_time=as.factor(lft_delivery_time),
         test_to_tracing=as.factor(test_to_tracing)) %>% 
  plotting_func(x=.,
                x_var=n_tests,
                group_var = lft_delivery_time,
                col_vars = test_to_tracing)

  a$plot+facet_nested(~test_to_tracing,labeller = labeller(lft_delivery_time=function(x)paste(x,"days postage delay"),test_to_tracing=function(x)paste(x,"days tracing delay")))+scale_color_viridis_d(option="magma",begin=0.2,end=0.8,name="LFT postage delay")+scale_y_continuous(limits=c(0.45,1),values = scales::percent)

  save_plot(dpi = 400, 
            device = "png",
            prefix = "delays",
            base = "_plot", 
            width = 210, 
            height = 100)  
  
#####################
#### RISK RATIOS ####
#####################

baseline <- get(results_name) %>%
  as_tibble() %>% 
  filter(quar_dur  == 10, 
         #test_to_tracing==1,
         is.na(assay)) %>%
  filter(!(is.infinite(inf_start) | is.infinite(inf_end))) %>% 
  select(ind_idx, max_overlap, inf_end, inf_start) %>% 
  group_by(ind_idx) %>% 
  summarise(all=sum(inf_end-inf_start),
            baseline_prop=sum(max_overlap)/all) %>% 
  mutate(name="baseline")


#### vs high baseline ----
get(results_name) %>% 
  mutate(test_to_tracing=as.factor(test_to_tracing)) %>% 
  rr_func(x=.,
          baseline = baseline,
          x_var = n_tests,
          col_vars = test_to_tracing,
          group_var = assay,
          log = T) 

fig3 <- a$plot+(b$plot+guides(colour=FALSE)+labs(y=""))

p_ranges_y <- c(10^(ggplot_build(fig3[[1]])$layout$panel_scales_y[[1]]$range$range),
                10^(ggplot_build(fig3[[2]])$layout$panel_scales_y[[1]]$range$range))

fig3+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")&ylim(min(p_ranges_y), max(p_ranges_y))

save_plot(dpi = 400, 
          device = "png",
          prefix = "baseline_high",
          base = "RR_plot", 
          width = 210, 
          height = 120)
