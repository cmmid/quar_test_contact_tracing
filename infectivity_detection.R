data.frame(x=-15:15) %>% 
  mutate(y=0.95*dgamma(x=x+infect_shift,shape=infect_shape,rate=infect_rate)/0.15) %>% 
  ggplot()+
  geom_line(aes(x=x,y=y,colour="Up-scaled infectivity profile * 0.95"))+
  geom_line(data=curves %>% 
              nest(-c(idx)) %>% 
              sample_n(100) %>% 
              unnest(),
            aes(x=days_since_infection-5,
                y=value,
                group=idx,
                colour="PCR curve posterior draws"),
            alpha=0.1)+
  labs(x="Days since onset",
       y="Assumed detection probability")+
  scale_colour_manual(name="",values = col_pal[2:3])+
  plotting_theme

save_plot(device="png",width=300,height=150)
