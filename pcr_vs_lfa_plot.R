results <- read_results(results_name)
results %>% filter(
  waning == "adhere_100",
  index_test_delay == 1,
  delay_scaling == 1,
  yvar == "infectivity_averted",
  quar_dur         %in% c(0,2,4,6,8,10,12,14),
  #type == "all"
  ) %>% 
  mutate(stringency=ifelse(stringency=="none","None",stringency))%>% 
  mutate(stringency=factor(stringency)) %>% 
  mutate(stringency=fct_relevel(stringency,"None","Tracing only", "End only", "Tracing and end")) %>% 
  mutate(assay=ifelse(stringency=="None","None",assay))%>% 
  mutate(assay=factor(assay)) %>% 
  mutate(assay=fct_relevel(assay,"None",after=0L)) %>% 
  ggplot(aes(x = time_since_exp, y = `50%`)) + 
  #RcmdrPlugin.KMggplot2::geom_stepribbon(aes(ymin = `2.5%`,
  #                                            ymax = `97.5%`,
  #                                            fill = assay),
  #                                        alpha = 0.5) +
  # geom_step(aes(color = assay),size=1) + 
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=assay),position=position_dodge(width=1),size=1,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=assay),position=position_dodge(width=1),size=1,alpha=0.5)+
  geom_point(aes(y = `50%`,
                 color = assay),
             #pch="-",
             size=1,
             position=position_dodge(width=1)) +
  scale_x_continuous(minor_breaks = breaks_width(2),
                     breaks       = breaks_width(2)
                     )+
  scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))+
  labs(x=expression("Quarantine required until"~italic("t")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_grid(#type 
    type ~ stringency,labeller = labeller(type=capitalize)) + 
  plotting_theme+
  scale_colour_manual(name="Assay",values = covid_pal)+
  scale_fill_manual(name="Assay",values = covid_pal)

save_plot(dpi = 400, 
          device = "png",
          prefix = "pcr_vs_lfa",
          base = "plot", 
          width = 210, 
          height = 140)

