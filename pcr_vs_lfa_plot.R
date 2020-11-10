results <- read_results(results_name)
results %>% filter(
  waning == "adhere_100",
  #index_test_delay == 1,
  #delay_scaling == 1,
  yvar == "infectivity_averted",
  quar_dur         %in% c(0,2,4,6,8,10,12,14),
  #type == "all"
  ) %>% 
  #mutate(stringency=fct_relabel(stringency,"Test upon tracing and end of quarantine","Test upon tracing\nand end of quarantine")) %>% 
  ggplot(aes(x = time_since_exp, y = `50%`)) + 
  #RcmdrPlugin.KMggplot2::geom_stepribbon(aes(ymin = `2.5%`,
  #                                            ymax = `97.5%`,
  #                                            fill = assay),
  #                                        alpha = 0.5) +
  # geom_step(aes(color = assay),size=1) + 
  geom_rect(data=. %>% filter(assay!="LFA"),aes(xmin=0,xmax=test_to_tracing*delay_scaling,ymin=0,ymax=1),alpha=0.05,fill="#A7A8AA")+
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
  #scale_y_continuous(limits = c(0,1),labels = scales::percent_format(accuracy = 1))+
  labs(x=expression("Quarantine required until"~italic("t")~"days have passed since exposure"),
       y="Transmission potential averted")+
  facet_nested(nest_line=T,
    delay_scaling+type ~ adherence+stringency,labeller = labeller(type=capitalize,delay_scaling=delay_scaling_labeller,
                                                    adherence=c("1"="100% self-isolate\nupon symptom onset",
                                                                "0.5"="50% self-isolate\nupon symptom onset",
                                                                "0"="0% self-isolate\nupon symptom onset"))) + 
  plotting_theme+
  scale_colour_brewer(name="Assay",palette = "Dark2")+
  scale_fill_brewer(name="Assay",palette="Dark2")

save_plot(dpi = 400, 
          device = "png",
          prefix = "pcr_vs_lfa_giant",
          base = "plot", 
          width = 500, 
          height = 400)

results %>% filter(
  waning == "adhere_100",
  #index_test_delay == 1,
  #delay_scaling == 1,
  yvar == "infectivity_averted",
  quar_dur         %in% c(0,2,4,6,8,10,12,14),
  #type == "all"
) %>% 
  #mutate(stringency=fct_relabel(stringency,"Test upon tracing and end of quarantine","Test upon tracing\nand end of quarantine")) %>% 
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
  facet_nested(nest_line=T,
               type ~ delay_scaling+stringency,labeller = labeller(type=capitalize,delay_scaling=delay_scaling_labeller)) + 
  plotting_theme+
  scale_colour_brewer(name="Assay",palette = "Dark2")+
  scale_fill_brewer(name="Assay",palette="Dark2")

save_plot(dpi = 400, 
          device = "png",
          prefix = "pcr_vs_lfa_all",
          base = "plot", 
          width = 300, 
          height = 200)

