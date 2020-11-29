baseline_scenario <- get(results_name) %>% 
  bind_rows() %>% 
  filter(adherence_iso==1,adherence_quar==1,delay_scaling==1,quar_dur==14,!tests) %>% 
  select(idx,sim,baseline_trans_pot_averted=trans_pot_averted) %>% 
  group_by(sim) %>% 
  summarise(baseline_prop=sum(baseline_trans_pot_averted)/n())

plot_1 <- all_results %>% 
  left_join(baseline_scenario) %>% 
  filter(#test_sensitivity==0.75,
  adherence_iso==1,
  adherence_quar==1,
  #delay_scaling==1,
  !multiple_tests
) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(sim,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(trans_pot_averted)/n()) %>% 
  left_join(baseline_scenario) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(quar_dur), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,
                     colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_log10(limits=c(0.4,NA))+
  labs(x=expression("Quarantine required until"~italic("n")~"days have passed since exposure"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed delays")+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme+
  scale_colour_manual(name="Strategy",values = col_pal[1:3])

plot_2 <- all_results %>% 
  left_join(baseline_scenario) %>% 
  filter(#test_sensitivity==0.75,
    adherence_iso==1,
    adherence_quar==1,
    #delay_scaling==1,
    multiple_tests
  ) %>%
  mutate(strategy=case_when(multiple_tests&tests~"Daily LFA testing",
                            tests&!multiple_tests&assay=="LFA"~"Post-exposure quarantine with LFA test",
                            tests&!multiple_tests&assay=="PCR"~"Post-exposure quarantine with PCR test",
                            !tests~"Post-exposure quarantine only"
  )) %>%
  group_by(sim,strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  summarise(prop=sum(trans_pot_averted)/n()) %>% 
  left_join(baseline_scenario) %>% 
  mutate(prop_ratio=prop/baseline_prop) %>% 
  group_by(strategy,adherence_quar,assay,quar_dur,n_tests,delay_scaling,sampling_freq) %>% 
  nest() %>%
  mutate(Q    = map(.x=data,
                    ~quantile(.$prop_ratio,
                              probs = c(0.025,0.25,0.5,0.75,0.975))),
         Mean = map_dbl(.x=data,
                        ~mean(.$prop_ratio,na.rm=T)),
         SD   = map_dbl(.x=data,
                        ~sd(.$prop_ratio,na.rm = T))) %>%
  unnest_wider(Q) %>% 
  mutate(strategy=factor(strategy)) %>% 
  ggplot(aes(x = factor(n_tests), y = `50%`)) + 
  geom_hline(aes(yintercept=1),linetype="dashed")+
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.3) +
  geom_linerange(aes(ymin = `25%`,
                     ymax = `75%`,colour=strategy),position=position_dodge(width=0.5),size=1.5,alpha=0.5)+
  geom_point(aes(y = `50%`,colour=strategy),
             #pch="-",
             size=1.5,
             position=position_dodge(width=0.5)) +
  # scale_x_continuous(#labels=delay_scaling_labeller,
  #                  guide=guide_axis(angle = 90))+
  # scale_x_continuous(minor_breaks = breaks_width(2),
  #                    breaks       = breaks_width(2)
  #)+
  scale_y_log10(limits=c(0.4,NA))+
  labs(x=expression("Daily LFA tests for"~italic("n")~"days after tracing"),
       y="Ratio of transmission potential averted compared to\nbaseline 14 day quarantine with observed delays")+
  scale_colour_manual(name="",values = col_pal[4])+
  facet_nested(nest_line=T,
               ~delay_scaling, labeller = labeller(
                 type = capitalize,
                 delay_scaling = delay_scaling_labeller,
                 adherence =
                   c("1" = "100% adhere\nto quarantine",
                     "0.5" =
                       "50% adhere\nto quarantine",
                     "0" =
                       "0% adhere\nto quarantine")
               )) +
  plotting_theme

plot_1+plot_2+plot_annotation(tag_levels = "A")+plot_layout(widths = c(3,2),guides = "collect")&theme(legend.position = "bottom")

save_plot(dpi = 600, 
          device = "png",
          prefix = "ratio",
          base = "plot", 
          width = 500, 
          height = 150)
