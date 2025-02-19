 peak_list<- unique(check_allother_peak$peak_key)


 for( i in 1:length(peak_list)){

 		current_peak = peak_list[i]

 		tmp_check =  check_traces %>% dplyr::filter(peak_key  == current_peak)
 		tmp_tau_decay = check %>% dplyr::filter(peak_key == current_peak)




 avgtrace_Plot<-ggplot(tmp_check, aes(x=normTime, y = dFF,  group=peak_key, colour=sensor))+
					                                                    #geom_vline(xintercept=0,lty='dashed',colour='black',size=0.2)+
					                                                    geom_line(alpha=0.8,size=1)+

					                                                #    geom_path(data = tmp_check, aes(x=normTime, y = .fitted, group=peak_key),colour='black',linetype="dashed", size=2,alpha=0.3)+
					                                                    
            
					                                                    
					                                                    labs( x="Normalized time (s)",
					                                                          y=expression(Delta*"F/F"),
					                                                          subtitle = paste0(current_peak, "___",tmp_tau_decay$tau_decay_ms)
					                                                          )+
					                                                    scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-HTL-AM")))+
															            scale_x_continuous(expand=c(0,0),breaks=c(0,0.5),oob=scales::oob_keep)+ #limits=c(-0.4,0.6)
															            scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1))+#,limits=c(-0.1,0.6))+
															            coord_cartesian(clip='off',xlim=c(-0.2,0.7),ylim =c(-0.2,1.0))+
															            
					                                                    theme_tufte()+
					                                                    my.theme+
					                                                    #facet_grid(sensor~chemical_condition, labeller = label_parsed,scale="free_y")+
					                                                    theme(strip.text.y=element_blank(),
					                                                    	  legend.position = "none",
					                                                    	  strip.text.x = element_blank())#,
					                                                    	  #plot.margin=unit(c(10,0,0,0), "cm"),
					                                                    	  #axis.line = element_blank(),
										                                      #axis.text.x = element_blank(),
										                                      #axis.text.y = element_blank(),
										                                      #axis.title = element_blank(),
										                                      #axis.ticks = element_blank())
avgtrace_Plot
					                                                    assign(paste0("check_tau_decay_avg_peak",i), avgtrace_Plot, envir = .GlobalEnv)

}