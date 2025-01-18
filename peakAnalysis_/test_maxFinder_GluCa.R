subDir <- "transmission_stats_revised_v9_ROIoutput"   


#dummy script to see if maxFinder_GluCa_v1.R can run within another function script. 

source(paste0(path,"sub-functions/setFigureDirectory.R")) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))
source(paste0(path,"peakAnalysis/maxFinder_GluCa.R"))
source(paste0(path,"peakAnalysis/peak_clipper.R"))
source(paste0(path,"peakAnalysis/png_plotter_v5_ROIoverlay_syntransmission.R"))
source(paste0(path,"peakAnalysis/png_plotter_v2_synTrans.R"))



library(ggforce)
library(scales)
library(ggthemes)
library(ggpmisc)
library(ggpubr)


transmission_stats<- function(traces_df = df_traces, saveROIs = ROIs_to_save){

						window_boundary = c(0.01,0.25)
						graph_output = FALSE
						plot_dim = 16
						scalefactor=0.75



						trace_times<- traces_df %>% group_by(chemical_condition, trackROI_key) %>% summarise(max_time = max(absoluteTime, na.rm=TRUE))
						trace_times$chemical_condition = factor(trace_times$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
						
						output_list<-maxFinder(traces_df, window_boundary, graph_output)

						max_y = 1.1
						max_x = 3
	
						amplitude_dist_full<- output_list[[1]]  
						amplitude_dist_full<- amplitude_dist_full %>% dplyr::filter(peak_color != "white")

						n_transmission <- amplitude_dist_full %>% dplyr::filter(peak_color != "baseline") %>%
																  group_by(chemical_condition, trackROI_key) %>%
																  summarise(n_total = n())
						frequency_dist <- left_join(n_transmission, trace_times) %>% 
												   group_by(chemical_condition, trackROI_key) %>% 
												   summarise(max_time = max_time,
												   			 #n_transmission = n(),
												   			  n_total = n_total,
												   			  freq_transmission = n_total/max_time) 

						my_avgs_full <- suppressMessages(
											amplitude_dist_full %>%
														group_by(chemical_condition,peak_color,trackROI_key) %>%
														mutate(event_interval = time_of_JFPeak - time_of_GluPeak,
																n_transmission  = n()) %>%
														ungroup() %>%
														group_by(chemical_condition, peak_color) %>%
														summarise(n = n(),
																	avg_n = mean(n_transmission,na.rm=TRUE),
																	sd_n = sd(n_transmission,na.rm=TRUE),


																  num_trackROI = length(unique(trackROI_key)),
																  
																  median_Glu_dFF = median(Glu_dFF, na.rm=TRUE),
																  mean_Glu_dFF = mean(Glu_dFF, na.rm=TRUE),
																  sd_Glu_dFF = sd(Glu_dFF, na.rm=TRUE),
																  se_Glu_dFF = sd_Glu_dFF/sqrt(n),

																  median_JF_dFF = median(JF_dFF, na.rm=TRUE),
																  mean_JF_dFF = mean(JF_dFF, na.rm=TRUE),
																  sd_JF_dFF = sd(JF_dFF, na.rm=TRUE),
																  se_JF_dFF = sd_JF_dFF/sqrt(n),

																  median_Glu_z_score =median(z_score_Glu, na.rm=TRUE), 
																  mean_Glu_z_score = mean(z_score_Glu, na.rm=TRUE),
																  sd_Glu_z_score = sd(z_score_Glu, na.rm=TRUE),
																  se_Glu_z_score = sd_Glu_z_score/sqrt(n),

																  median_JF_z_score =median(z_score_JF, na.rm=TRUE), 
																  mean_JF_z_score = mean(z_score_JF, na.rm=TRUE),
																  sd_JF_z_score = sd(z_score_JF, na.rm=TRUE),
																  se_JF_z_score = sd_JF_z_score/sqrt(n),


																  median_event_interval = median(event_interval, na.rm=TRUE),
																  mean_event_interval = mean(event_interval, na.rm=TRUE),
																  sd_event_interval = sd(event_interval, na.rm=TRUE),
																  se_event_interval = sd_event_interval/sqrt(n)
																  ) 
											)


							ROI_avgs <- suppressMessages(
											amplitude_dist_full %>%
														group_by(chemical_condition,peak_color,trackROI_key) %>%
														mutate(event_interval = time_of_JFPeak - time_of_GluPeak,
																n_transmission  = n()) %>%
														ungroup() %>%
														group_by(chemical_condition, trackROI_key) %>%
														dplyr::filter(peak_color == "blue") %>%
														summarise(n = n(),
																	avg_n = unique(n_transmission),
																	#sd_n = #sd(n_transmission,na.rm=TRUE),


																  #num_trackROI = length(unique(trackROI_key)),
																  corr_slope = coef(lm(z_score_JF ~ z_score_Glu))[[2]],
																  median_Glu_dFF = median(Glu_dFF, na.rm=TRUE),
																  mean_Glu_dFF = mean(Glu_dFF, na.rm=TRUE),
																  sd_Glu_dFF = sd(Glu_dFF, na.rm=TRUE),
																  se_Glu_dFF = sd_Glu_dFF/sqrt(n),

																  var_Glu_dFF = var(Glu_dFF,na.rm=TRUE),

																  median_JF_dFF = median(JF_dFF, na.rm=TRUE),
																  mean_JF_dFF = mean(JF_dFF, na.rm=TRUE),
																  sd_JF_dFF = sd(JF_dFF, na.rm=TRUE),
																  se_JF_dFF = sd_JF_dFF/sqrt(n),

																  var_JF_dFF = var(JF_dFF, na.rm=TRUE),

																  median_Glu_z_score =median(z_score_Glu, na.rm=TRUE), 
																  mean_Glu_z_score = mean(z_score_Glu, na.rm=TRUE),
																  sd_Glu_z_score = sd(z_score_Glu, na.rm=TRUE),
																  se_Glu_z_score = sd_Glu_z_score/sqrt(n),

																  var_Glu_z_score = var(z_score_Glu, na.rm=TRUE),

																  median_JF_z_score =median(z_score_JF, na.rm=TRUE), 
																  mean_JF_z_score = mean(z_score_JF, na.rm=TRUE),
																  sd_JF_z_score = sd(z_score_JF, na.rm=TRUE),
																  se_JF_z_score = sd_JF_z_score/sqrt(n),

																  var_JF_z_score = var(z_score_JF, na.rm=TRUE),

																  median_event_interval = median(event_interval, na.rm=TRUE),
																  mean_event_interval = mean(event_interval, na.rm=TRUE),
																  sd_event_interval = sd(event_interval, na.rm=TRUE),
																  se_event_interval = sd_event_interval/sqrt(n)
																  ) %>%
														dplyr::filter(mean_event_interval < 0.5,
																		corr_slope < 10,
																		corr_slope > -10 ) 
											)

						

						my_medians_full<- suppressMessages(
											amplitude_dist_full %>%
														group_by(chemical_condition,peak_color) %>%
														#dplyr::filter(peak_color != "white") %>% 
														summarise(median_JF_dFF = median(JF_dFF, na.rm=TRUE),
																  median_Glu_dFF = median(Glu_dFF, na.rm=TRUE),

																  median_JF_z_score = median(z_score_JF, na.rm=TRUE),
																  median_Glu_z_score = median(z_score_Glu, na.rm=TRUE)
																  ) 
											)




						#levels=c("ctrl","vehicle","APV")
						amplitude_dist_slim<- amplitude_dist_full %>% ungroup() %>% 
																	 dplyr::filter(Glu_dFF <= max_x, JF_dFF <= max_y, peak_color == "blue") %>%
																	 group_by(chemical_condition,peak_color) %>%
																	 mutate(peak_color_revised = case_when(peak_color == "baseline" ~ "baseline",
																	 										peak_color == "blue" & chemical_condition == "Pre-Treatment" ~ "blue",
																	 										peak_color == "blue" & chemical_condition == expression("100"~mu*"M AP5") ~ "black",
																	 										peak_color == "blue" & chemical_condition == "Control" ~ "lightblue") 
																	 										) #%>%
																	 #dplyr::filter()

						my_medians_slim<- suppressMessages(
											amplitude_dist_slim %>%
														group_by(chemical_condition,peak_color) %>%
														summarise(median_JF_dFF = median(JF_dFF, na.rm=TRUE),
																  median_Glu_dFF = median(Glu_dFF, na.rm=TRUE)
																  ) 
											)

						print("Attempting amplitude_scatterplot on the full dataset.")

						 formula = y ~ x

						amplitude_scatterplot<-ggplot(amplitude_dist_slim, aes(x=Glu_dFF, y=JF_dFF, group=chemical_condition,colour=peak_color))+#,alpha=peak_color))+
		                                                geom_point(size=4,alpha=0.6)+
		                                                geom_smooth(data=subset(amplitude_dist_full,peak_color == "blue"), method = "lm",formula = y~ x, se = FALSE, colour = "black", size=1.5,alpha=0.8)+
		                                                stat_cor( method = "pearson",color = "black", geom = "text",label.x = 0.3*max_x,label.y=0.95*max_y,size=7)+ 
		                                                stat_regline_equation(label.x = 0.3*max_x, label.y = 0.85*max_y,color="black",size=7) +
		                                                labs( x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
		                                                      y=expression(JF[646]*" amplitude ("*Delta*"F/F)"),
		                                                  	  subtitle="Spine Population Behavior",
		                                                  	  tag="H")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
		                                                scale_colour_manual(name = "", breaks=c("blue","white", "baseline"),
		                                           									values = c("blue"='blue',"white" = "grey31","baseline" = "grey61"),
		                                           									labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"))+
		                                                #scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
		                                                coord_cartesian(xlim=c(0,max_x),ylim=c(-0.1,max_y))+
		                                                scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
		                                                scale_x_continuous(breaks=c(0,1,2,3))+
		                                                facet_grid(~chemical_condition,labeller = label_parsed)+
		                                                theme_tufte()+
		                                                my.theme +
		                                                theme(legend.position = "right")#,
		                                           		      #strip.text.x = strip_text_arg_x)
		                save_plot(filename=paste0("amplitude_corr_full_dataset"),plot=amplitude_scatterplot, plot_width=plot_dim*2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)


		                max_x_z_score = 5#max(amplitude_dist_slim$z_score_Glu, na.rm=TRUE)*1.2
		                max_y_z_score = 5#max(amplitude_dist_slim$z_score_JF, na.rm=TRUE)*1.2

		                min_x_z_score = -1.5#min(amplitude_dist_slim$z_score_Glu, na.rm=TRUE)*1.2
		                min_y_z_score = -1.5#min(amplitude_dist_slim$z_score_JF, na.rm=TRUE)*1.2

		                break_vec_x = c(0,2,4)
		                break_vec_y = c(0,2,4)
		                
		                amplitude_dist_slim_z_score <- amplitude_dist_slim %>% dplyr::filter(z_score_Glu <= max_x_z_score, z_score_Glu >= min_x_z_score, z_score_JF <= max_y_z_score, z_score_JF >= min_y_z_score)

		       


		                amplitude_z_score_scatterplot<-ggplot(amplitude_dist_slim_z_score, aes(x=z_score_Glu, y=z_score_JF, group=chemical_condition,colour=peak_color_revised))+#,alpha=peak_color))+
		                                                #geom_hline(yintercept = c(0.246,0.0311,-0.461))+
		                                                geom_point(size=4,alpha=0.45)+
		                                                stat_smooth(colour="black",method = "lm",formula = formula, se = TRUE, size=1.5,alpha=1)+
		                                                stat_regline_equation(formula = formula,  label.x = min_x_z_score*0.5, label.y = 0.82*max_y_z_score,color="black",size=7) +
		                                                stat_cor( method = 'pearson',color = "black", geom = "text",label.x = min_x_z_score*0.5,label.y=0.95*max_y_z_score,size=7,p.accuracy = 0.001, r.accuracy = 0.01)+ 
		                                                
		                                                labs( x=expression("iGluSnFR3 signal z-score"),
		                                                      y=expression(JF[646]*" signal z-score"),
		                                                  	  #subtitle="Population Dose-Response",
		                                                  	  tag="E")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
		                                                scale_colour_manual(name = "", breaks=c("blue","black","lightblue", "baseline"),
		                                           										values = c("blue"='blue',"black"="grey31","lightblue" = "skyblue3","baseline" = "grey61"),
		                                           									labels=c("blue" = "Transmission", "black" = "APV","lightblue" = "Control", "baseline" = "Baseline"))+
		                                           									#labels=c("blue" = "Pre-Treatment", "black" = "APV","lightblue" = "Control","baseline" = "Baseline"))+
		                                                #scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
		                                                coord_cartesian(xlim=c(min_x_z_score,max_x_z_score),ylim=c(min_y_z_score,max_y_z_score))+
		                                               	scale_y_continuous(breaks = break_vec_y)+
		                                               	scale_x_continuous(breaks = break_vec_x)+
		                                                facet_grid(~chemical_condition,labeller = label_parsed)+
		                                                theme_tufte()+
		                                                my.theme +
		                                                theme(legend.position = "none")#,
		                                           		      #strip.text.x = strip_text_arg_x)
		                save_plot(filename=paste0("amplitude_corr_z_score_full_dataset"),plot=amplitude_z_score_scatterplot, plot_width=plot_dim*2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)

						amplitude_z_score_corr<<- amplitude_z_score_scatterplot



### mean vs variance? 


		                amplitude_z_score_scatterplot<-ggplot(ROI_avgs, aes(x=mean_JF_z_score, y=var_JF_z_score, group=chemical_condition))+#,alpha=peak_color))+
		                                                #geom_hline(yintercept = c(0.246,0.0311,-0.461))+
		                                                geom_point(size=4,alpha=0.45)+
		                                                stat_smooth(colour="black",method = "lm",formula = formula, se = TRUE, size=1.5,alpha=1)+
		                                                stat_regline_equation(formula = formula,  label.x = min_x_z_score*0.5, label.y = 0.85*max_y_z_score,color="black",size=7) +
		                                                stat_cor( method = 'pearson',color = "black", geom = "text",label.x = min_x_z_score*0.5,label.y=0.95*max_y_z_score,size=7,p.accuracy = 0.001, r.accuracy = 0.01)+ 
		                                                
		                                                labs( x=expression(JF[646]*" z-score"),
		                                                      y=expression("Variance of "*JF[646]*" z-score"),
		                                                  	  #subtitle="Population Dose-Response",
		                                                  	  tag="E")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
		                                                #scale_colour_manual(name = "", breaks=c("blue","black","lightblue", "baseline"),
		                                           	#									values = c("blue"='blue',"black"="grey31","lightblue" = "skyblue3","baseline" = "grey61"),
		                                           #									labels=c("blue" = "Transmission", "black" = "APV","lightblue" = "Control", "baseline" = "Baseline"))+
		                                           									#labels=c("blue" = "Pre-Treatment", "black" = "APV","lightblue" = "Control","baseline" = "Baseline"))+
		                                                #scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
		                                                #coord_cartesian(xlim=c(min_x_z_score,max_x_z_score),ylim=c(min_y_z_score,max_y_z_score))+
		                                               	#scale_y_continuous(breaks = break_vec_y)+
		                                               	#scale_x_continuous(breaks = break_vec_x)+
		                                                facet_grid(~chemical_condition,labeller = label_parsed)+
		                                                theme_tufte()+
		                                                my.theme +
		                                                theme(legend.position = "none")#,
		                                           		      #strip.text.x = strip_text_arg_x)
		                save_plot(filename=paste0("amplitude_var_correlation_z_score_full_dataset"),plot=amplitude_z_score_scatterplot, plot_width=plot_dim*2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)



		                high_var_ROIs<- ROI_avgs %>% dplyr::filter(chemical_condition == "Pre-Treatment", var_JF_z_score >= quantile(var_JF_z_score, 0.5, na.rm=TRUE))
		                print(paste0("The total number of high variance spines is: ", length(unique(high_var_ROIs$trackROI_key))))

		                subset_amplitude_dist_slim_z_score <- amplitude_dist_slim_z_score %>% dplyr::filter(trackROI_key %in% unique(high_var_ROIs$trackROI_key))


		                amplitude_z_score_scatterplot<-ggplot(subset_amplitude_dist_slim_z_score, aes(x=z_score_Glu, y=z_score_JF, group=chemical_condition,colour=peak_color_revised))+#,alpha=peak_color))+
		                                                #geom_hline(yintercept = c(0.246,0.0311,-0.461))+
		                                                geom_point(size=4,alpha=0.45)+
		                                                stat_smooth(colour="black",method = "lm",formula = formula, se = TRUE, size=1.5,alpha=1)+
		                                                stat_regline_equation(formula = formula,  label.x = min_x_z_score*0.5, label.y = 0.85*max_y_z_score,color="black",size=7) +
		                                                stat_cor( method = 'pearson',color = "black", geom = "text",label.x = min_x_z_score*0.5,label.y=0.95*max_y_z_score,size=7,p.accuracy = 0.001, r.accuracy = 0.01)+ 
		                                                
		                                                labs( x=expression("iGluSnFR3 signal z-score"),
		                                                      y=expression(JF[646]*" signal z-score"),
		                                                  	  #subtitle="Population Dose-Response",
		                                                  	  tag="E")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
		                                                scale_colour_manual(name = "", breaks=c("blue","black","lightblue", "baseline"),
		                                           										values = c("blue"='blue',"black"="grey31","lightblue" = "skyblue3","baseline" = "grey61"),
		                                           									labels=c("blue" = "Transmission", "black" = "APV","lightblue" = "Control", "baseline" = "Baseline"))+
		                                           									#labels=c("blue" = "Pre-Treatment", "black" = "APV","lightblue" = "Control","baseline" = "Baseline"))+
		                                                #scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
		                                                coord_cartesian(xlim=c(min_x_z_score,max_x_z_score),ylim=c(min_y_z_score,max_y_z_score))+
		                                               	scale_y_continuous(breaks = break_vec_y)+
		                                               	scale_x_continuous(breaks = break_vec_x)+
		                                                facet_grid(~chemical_condition,labeller = label_parsed)+
		                                                theme_tufte()+
		                                                my.theme +
		                                                theme(legend.position = "none")#,
		                                           		      #strip.text.x = strip_text_arg_x)
		                save_plot(filename=paste0("amplitude_corr_z_score_high_Var_ROIs"),plot=amplitude_z_score_scatterplot, plot_width=plot_dim*2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)

						#amplitude_z_score_corr<<- amplitude_z_score_scatterplot


						high_corr_ROIs<- ROI_avgs %>% dplyr::filter(chemical_condition == "Pre-Treatment", corr_slope >= quantile(corr_slope, 0.5, na.rm=TRUE))
		                


						 amplitude_z_score_scatterplot<-ggplot(subset_amplitude_dist_slim_z_score, aes(x=z_score_Glu, y=z_score_JF, group=chemical_condition,colour=peak_color_revised))+#,alpha=peak_color))+
		                                                #geom_hline(yintercept = c(0.246,0.0311,-0.461))+
		                                                geom_point(size=4,alpha=0.45)+
		                                                stat_smooth(colour="black",method = "lm",formula = formula, se = TRUE, size=1.5,alpha=1)+
		                                                stat_regline_equation(formula = formula,  label.x = min_x_z_score*0.5, label.y = 0.85*max_y_z_score,color="black",size=7) +
		                                                stat_cor( method = 'pearson',color = "black", geom = "text",label.x = min_x_z_score*0.5,label.y=0.95*max_y_z_score,size=7,p.accuracy = 0.001, r.accuracy = 0.01)+ 
		                                                
		                                                labs( x=expression("iGluSnFR3 signal z-score"),
		                                                      y=expression(JF[646]*" signal z-score"),
		                                                  	  #subtitle="Population Dose-Response",
		                                                  	  tag="E")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
		                                                scale_colour_manual(name = "", breaks=c("blue","black","lightblue", "baseline"),
		                                           										values = c("blue"='blue',"black"="grey31","lightblue" = "skyblue3","baseline" = "grey61"),
		                                           									labels=c("blue" = "Transmission", "black" = "APV","lightblue" = "Control", "baseline" = "Baseline"))+
		                                           									#labels=c("blue" = "Pre-Treatment", "black" = "APV","lightblue" = "Control","baseline" = "Baseline"))+
		                                                #scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
		                                                coord_cartesian(xlim=c(min_x_z_score,max_x_z_score),ylim=c(min_y_z_score,max_y_z_score))+
		                                               	scale_y_continuous(breaks = break_vec_y)+
		                                               	scale_x_continuous(breaks = break_vec_x)+
		                                                facet_grid(~chemical_condition,labeller = label_parsed)+
		                                                theme_tufte()+
		                                                my.theme +
		                                                theme(legend.position = "none")#,
		                                           		      #strip.text.x = strip_text_arg_x)
		                save_plot(filename=paste0("amplitude_corr_z_score_high_correlation_ROIs"),plot=amplitude_z_score_scatterplot, plot_width=plot_dim*2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)

						#amplitude_z_score_corr<<- amplitude_z_score_scatterplot



						###Determine whether transmission stats as calculated here yield an interesting frequency plot - the answer is no at the moment! 

						bar_data<- amplitude_dist_full %>% group_by(chemical_condition, peak_color) %>% 
														   dplyr::filter(peak_color != "baseline") %>% 
														   summarise(n=n()) %>% 
														   mutate(freq=n/sum(n)) %>% 
														   mutate(asPercent = label_percent()(freq),
						                                          groupfact = "All Spines")


						#for printing
						check_bar_data<-bar_data %>% select(chemical_condition,peak_color,n,asPercent)
						print(check_bar_data)
						
						

						
						amplitude_dist_full$peak_color = factor(amplitude_dist_full$peak_color,levels=c("baseline","white","blue"))
						amplitude_dist_slim$peak_color = factor(amplitude_dist_slim$peak_color,levels=c("baseline","white","blue"))
											

						# Function to create histogram plot
						create_histogram_plot <- function(data,median_data, x_variable, median_var, x_label, filename, lower_bound_x, upper_bound_x,color_override,global_output,tag_var,legend_bool) {

								if(upper_bound_x == 0.8){
								  	break_vec = c(0,0.4,0.8)
								  } else if(upper_bound_x == 3){
								  	break_vec = c(0,1,2,3)
								  } else {
								  	break_vec = seq(0,upper_bound_x,length.out=5)
								  }
								  print(paste0("creating a histogram for: ", x_variable) )
								  bin = (upper_bound_x - lower_bound_x) / 50
								 
								  plot <- ggplot(data, aes(x = get(x_variable), group = peak_color, fill = peak_color, colour = peak_color)) +
								    geom_vline(data = median_data, aes(xintercept = get(median_var), colour = peak_color), size = 2, lty = "solid") +
								    geom_freqpoly(aes(y = ..ncount.., colour=peak_color), binwidth = bin, alpha = 0.7, boundary = 0, size = 1) +
								    geom_histogram(aes(y = ..ncount.., fill=peak_color), binwidth = bin, alpha = 0.5, boundary = 0, position = "identity", colour = "black", size = 0.8) +
								    {if(legend_bool)geom_text(data=median_data, aes(label = chemical_condition), x = upper_bound_x+(0.05*upper_bound_x), y = 0.45, colour="black",size=8,hjust=1,check_overlap = TRUE,parse=TRUE,show_guide=FALSE)}+
		    
								    labs(
								      x = x_label,
								      y = expression(N/N[max]), #expression("Density"),#
								      fill = "",
								      tag = tag_var
								    ) +
								    scale_colour_manual(name = "", breaks=c("blue","white", "baseline"),
				                                           		   values = c("blue"=color_override,"white" = "grey31","baseline" = "grey61"),
				                                           		   labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"),
				                                           		   guide="none")+
								    
								    scale_fill_manual(name = "", breaks=c("blue","white", "baseline"),
		 														 values = c("blue"=color_override,"white" = "grey31","baseline" = "grey61"),
				                                           		 labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"))+
								    coord_cartesian(ylim = c(0, 1.3), xlim = c(lower_bound_x, upper_bound_x)) +
								    scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
								    scale_x_continuous(breaks=break_vec)+
								    theme_tufte() +
								    facet_grid(chemical_condition ~ ., labeller = label_parsed) +
								    my.theme +
								    theme(
								      legend.position = "top",
								      legend.direction = "vertical",
								      legend.justification = "center",
								      legend.title = element_blank(),
								      #axis.text.y = axis_text_y_arg,
								      #axis.title.y = axis_title_y_arg,
								      strip.text.y.right = element_blank()#element_text(colour="black", size=scalefactor*34, family="sans", angle=0)
								    )

								  save_plot(
								    filename = paste0(filename, ""),
								    plot = plot,
								    plot_width=plot_dim*1.5,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

								  if(global_output == TRUE){
								  	assign(paste0(x_variable, "_histo"), plot, envir = .GlobalEnv)
								  }
						}

					
						create_histogram_plot(amplitude_dist_full, my_medians_full, "Glu_dFF", "median_Glu_dFF", expression("iGluSnFR3 amplitude ("*Delta*"F/F)"), "histo___amplitude_dist_full____Glu_amplitude", -0.15, 2,"forestgreen",TRUE,"H",TRUE)
						create_histogram_plot(amplitude_dist_full, my_medians_full, "JF_dFF", "median_JF_dFF", expression(JF[646]*" amplitude ("*Delta*"F/F)"), "histo___amplitude_dist_full____JF646_amplitude", -0.15, 0.8,"#FF33FF",TRUE,"I",FALSE)
						create_histogram_plot(amplitude_dist_full, my_medians_full, "z_score_JF", "median_JF_z_score", expression(JF[646]*" signal z-score"), "histo___amplitude_dist_full____JF646_z-score", -1.5, 5,"#FF33FF",TRUE,"I",FALSE)
						create_histogram_plot(amplitude_dist_full, my_medians_full, "z_score_Glu", "median_Glu_z_score", expression("iGluSnFR3 signal z-score"), "histo___amplitude_dist_full____Glu_z-score", -1.5, 5,"forestgreen",TRUE,"I",FALSE)



							
						create_violin_plot<- function(data, y_variable,  y_label, colour_arg, tag_var,  filename){							  

								 print(paste0("creating a violin plot for: ", y_variable) )
								 

							  plot <- ggplot(data, aes(x = chemical_condition, y = get(y_variable))) +
							    geom_sina(colour = colour_arg,size=1.5,alpha=0.5)+
							    geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
							    

							    labs(
							      y = y_label,
							      x = "",
							      tag = tag_var
							    ) +
							    theme_tufte() +
							    my.theme +
							    theme(legend.position = "none")

							  save_plot(filename = paste0(filename, ""), plot = plot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)
							}



ROI_avgs_pretreat <- ROI_avgs %>% dplyr::filter(chemical_condition == "Pre-Treatment")

create_violin_plot(ROI_avgs, "var_JF_z_score", expression("Variance of "*JF[646]~"z-score"), "#FF33FF", "X", "violin__var_JF_z_score")
create_violin_plot(ROI_avgs, "var_Glu_z_score", expression("Variance of iGluSnFR3 z-score"), "forestgreen", "X", "violin__var_Glu_z_score")
create_violin_plot(ROI_avgs, "corr_slope", expression("Slope of amplitude correlation"), "blue", "X", "violin__amplitude_correlation")
create_violin_plot(ROI_avgs, "mean_event_interval", expression(Delta*"t between pre- and postsynaptic events (s)"), "blue", "X", "violin__deltaT")













						#For stats testing
						amplitude_statistics<-amplitude_dist_full  %>% dplyr::filter(peak_color == "blue")
						print(amplitude_statistics)
						#AtoB<- amplitude_statistics %>% dplyr::filter(chemical_condition == "Pre-Treatment")
						#AtoC
						B<- amplitude_statistics %>% dplyr::filter(chemical_condition == "Control")
						C<- amplitude_statistics %>% dplyr::filter(chemical_condition == "\"100\" ~ mu * \"M AP5\"")

						AtoB<- amplitude_statistics %>% dplyr::filter(chemical_condition == "Pre-Treatment", trackROI_key %in% unique(B$trackROI_key))
						AtoC<- amplitude_statistics %>% dplyr::filter(chemical_condition == "Pre-Treatment", trackROI_key %in% unique(C$trackROI_key))
						
						

						ks_test<- wilcox.test(AtoB$Glu_dFF,B$Glu_dFF)
                        print("testing Glu_dFF, ctrl vs vehicle with wilcox.test")
                        print(ks_test) 
                         
                        ks_test<- wilcox.test(AtoC$Glu_dFF,C$Glu_dFF)
                        print("testing Glu_dFF, ctrl vs APV with wilcox.test")
                        print(ks_test) 
                        
                         ks_test<- wilcox.test(B$Glu_dFF,C$Glu_dFF)
                        print("testing Glu_dFF, vehicle vs APV with wilcox.test")
                        print(ks_test) 
                        
                        

						ks_test<- wilcox.test(AtoB$JF_dFF,B$JF_dFF)
                        print("testing JF_dFF, ctrl vs vehicle with wilcox.test")
                        print(ks_test) 
                         
                        ks_test<- wilcox.test(AtoC$JF_dFF,C$JF_dFF)
                        print("testing JF_dFF, ctrl vs APV with wilcox.test")
                        print(ks_test) 
                        
                         ks_test<- wilcox.test(B$JF_dFF,C$JF_dFF)
                        print("testing JF_dFF, vehicle vs APV with wilcox.test")
                        print(ks_test) 
                        


                        
						window_boundary = c(0.1,0.6)
						traces_binned<-peak_clipper(traces_df, window_boundary)
						print(unique(traces_binned$sensor))
						interFrame.var = as.vector(traces_df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE) 
						traces_binned<- traces_binned %>% group_by(sensor,chemical_condition, trackROI_key, peak_key) %>% dplyr::filter(max_dFF <= 3)
                        min_val = -0.05


                        traces_binned$chemical_condition = factor(traces_binned$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															
                         custom_y <- tribble(~sensor, ~absoluteTime, ~dFF,
											  "GluSnFR3", 0, min_val,
											  "GluSnFR3", 0, 0.4,
											  "JF646", 0, min_val,
											  "JF646", 0.1, 0.2
                  								)

                         #GluSnFR3 x segment,
                         #JF646 x_segment,
                         
                         #GluSnFR3 y segment,
                         #JF646 y_segment,
                         x_origin = -0.15 
                         x_final = 0.2+x_origin
                         x_length = x_final - x_origin

                         y_origin = 0.025
                         y_final = 0.2 + y_origin
                         y_length = y_final-y_origin

                         scalebar_x_data<-data.frame(sensor = c("GluSnFR3","JF646"),
                         						  chemical_condition = c('ctrl','ctrl'),
                         							x_start = c(NA,x_origin),
                         							x_end = c(NA,x_final),
                         							y_start = c(NA,y_origin),
                         							y_end = c(NA,y_origin)
                         							#x_expr = c(NULL,"200 ms"))#,
                         							#y_expr = c(expression("0.2 "*Delta[F/F]),expression("0.2 "*Delta[F/F]))
                         							)
                         x_label<- data.frame(sensor = c('JF646'),
                         						chemical_condition = c('ctrl'),
                         						x= (x_length/2)+x_origin,
                         						y=y_origin - 0.03)

                         scalebar_y_data<-data.frame(sensor = c("GluSnFR3","JF646"),
                         						  chemical_condition = c('ctrl','ctrl'),
                         							x_start = c(x_origin,x_origin),
                         							x_end = c(x_origin,x_origin),
                         							y_start = c(y_origin,y_origin),
                         							y_end = c(y_final,y_final)#,
                         							#y_expr = c(expression("0.2 "*Delta[F/F]),expression("0.2 "*Delta[F/F]))
                         							)

                         y_label <-data.frame(sensor = c('GluSnFR3','JF646'),
                         						chemical_condition = c('ctrl','ctrl'),
                         						x= c(x_origin - 0.2, x_origin - 0.2),
                         						y=c( (y_length/2)+y_origin,(y_length/2)+y_origin))


						 my_y_breaks <- function(x) {if(max(x) <= 0.4)  seq(0.0,0.3,by=0.1) else if(max(x) <= 1.5 & max(x) >= 0.5) seq(0.0,0.6,by=0.3) else seq(0.0,7.5,by=2.5)}
						 #my_breaks <- function(x) { if(max(x) < 20)  seq(0,15,5) else seq(0,25,5) }
						 formatter <- function(...){
													  function(x) format(round(x, 1), ...)
													}


						 traces_binned$sensor = as.factor(as.character(traces_binned$sensor))
						 levels(traces_binned$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))
					

						 custom_y$sensor = as.factor(as.character(custom_y$sensor))
						 levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))

						 scalebar_x_data$chemical_condition = factor(scalebar_x_data$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
						 scalebar_x_data$sensor = as.factor(as.character(scalebar_x_data$sensor))
						 levels(scalebar_x_data$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))
	 
						 scalebar_y_data$chemical_condition = factor(scalebar_y_data$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
						 scalebar_y_data$sensor = as.factor(as.character(scalebar_y_data$sensor))
						 levels(scalebar_y_data$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))

						 x_label$chemical_condition = factor(x_label$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
						 x_label$sensor = as.factor(as.character(x_label$sensor))
						 levels(x_label$sensor) <- c(expression(JF[646]*"-BAPTA-HTL-AM"))

						 y_label$chemical_condition = factor(y_label$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
						 y_label$sensor = as.factor(as.character(y_label$sensor))
						 levels(y_label$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))
	 
	 					
						 traces_binned <- traces_binned %>% dplyr::filter(normTime >= -0.1, normTime <= 0.6)#, peak_color == "blue")

						 count_peak_keys = length(unique(traces_binned$peak_key))
						 print(paste0("The total number of unique peaks that went into this visualization was : ", count_peak_keys ))
						 check_numbers<- my_avgs_full %>% ungroup() %>% group_by(chemical_condition) %>% dplyr::filter(peak_color != "baseline")
						 sum_events <- sum(check_numbers[ , "n"], na.rm=TRUE)
						 print(paste0("The total number of unique transmission events that went into the average calculations was : ", sum_events))

						         text_size = 9.09

					    avgtrace_Plot<-ggplot(traces_binned, aes(x=normTime, y = dFF,  group=peak_key, colour=sensor))+
					                                                    #geom_vline(xintercept=0,lty='dashed',colour='black',size=0.2)+
					                                                    stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), geom="line", fun=mean, binwidth=interFrame.var,size=0.5,alpha=1)+
					                                                    stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), geom = "errorbar", fun.data=mean_se, binwidth=interFrame.var,size=0.5,alpha=0.65)+
						                            
					                                                    stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), geom="point", fun=mean, binwidth=interFrame.var,size=0.8,alpha=1)+
					                                                    
					                                                    
					                                                    #stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), colour="grey51",geom="smooth", binwidth=interFrame.var,size=1,alpha=0.6)+
					                                                    #geom_path(size=0.4,alpha=0.1,colour="grey51")+
					                                                    geom_blank(data=custom_y,aes(absoluteTime,dFF), inherit.aes=FALSE)+
					                                                    geom_segment(data=scalebar_x_data, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 2.25, colour="black", inherit.aes=FALSE)+
															            geom_segment(data=scalebar_y_data, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 2.25, colour="black", inherit.aes=FALSE)+
															            geom_label(data=x_label,aes(x=x,y=y), label=expression('200 ms'), size= text_size, fill=NA,colour="black",label.size=NA,family="sans",inherit.aes=FALSE)+ #
 															            geom_label(data=y_label,aes(x=x,y=y), label=expression('0.2 '*Delta*"F/F"), size= text_size, fill=NA,colour="black",label.size=NA,family="sans",inherit.aes=FALSE)+ #
 															            
 															            #{if(add_scalebars)geom_label(aes(x=label_y$x,y=label_y$y, label=label_y_expr), nudge_x = (scalebar_x_length*-0.7), size= 7, colour="black",parse=TRUE,label.size=NA,family="sans")}+
															            
           
															            
					                                                    
					                                                    labs( x="Normalized time (s)",
					                                                          y=expression(Delta*"F/F"),
					                                                          tag = "D",
					                                                          subtitle = "Averaged Transmission Events"
					                                                          )+
					                                                    scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-HTL-AM")))+
															            scale_x_continuous(expand=c(0,0),breaks=c(0,0.5),oob=scales::oob_keep)+ #limits=c(-0.4,0.6)
															            scale_y_continuous(expand=c(0,0))+#,limits=c(-0.1,0.6))+
															            coord_cartesian(clip='off',xlim=c(-0.4,0.6))+
															            
					                                                    theme_tufte()+
					                                                    my.theme+
					                                                    facet_grid(sensor~chemical_condition, labeller = label_parsed,scale="free_y")+
					                                                    theme(strip.text.y=element_blank(),
					                                                    	  legend.position = "none",
					                                                    	  #plot.margin=unit(c(10,0,0,0), "cm"),
					                                                    	  axis.line = element_blank(),
										                                      axis.text.x = element_blank(),
										                                      axis.text.y = element_blank(),
										                                      axis.title = element_blank(),
										                                      axis.ticks = element_blank())

					                                                    assign("avg_traces_dualcolor", avgtrace_Plot, envir = .GlobalEnv)




						save_plot(filename=paste0("averagePeaks_dualColor_clipped"),plot=avgtrace_Plot, plot_width=plot_dim*1.2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)
					                         





## plotting avg waveforms by ROI


allROIs <- unique(traces_binned$trackROI_key)
n_transmission > 5




library(devtools)
library(nls.multstart) 
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)

drops = c('data', '.resid','fit',"(weights)") #.fitted
					


### try decay fitting routine on dataset
expDecayFit<- function(df) {			

					expFit<- df %>%
							dplyr::filter(sensor != "iGluSnFR3", chemical_condition == "Pre-Treatment") %>% #unit test ###, peak_key != "ctrl-dish10-plate04-region02-ROI0004-4", peak_key != "ctrl-dish10-plate04-region02-ROI0010-9", peak_key != "ctrl-dish10-plate04-region02-ROI0012-5"
							group_by(sensor,chemical_condition,trackROI_key,peak_key) %>%
							dplyr::filter(!is.na(dFF)) %>%
							slice(which.max(dFF):length(dFF)) %>%
							mutate(count_obs = n()) %>%
							dplyr::filter(count_obs >= 3) %>%
							#mutate(decay_normTime = normTime - min(normTime))  %>%
							#dplyr::filter(normTime < 1.0) %>%
		    					ungroup() %>%
		    					group_by(sensor,chemical_condition,trackROI_key,peak_key) %>%
								nest() %>%
		   					  	mutate(
		   						        fit = purrr::map(data, ~nls_multstart(dFF ~ A*exp(-normTime/tau_decay),
		   						        				data=.x,
		   						        				iter = 1000,
		   						        				start_lower = c(A = 0, tau_decay= 0.01 ),
		   						        				start_upper = c(A = 6, tau_decay= 0.7),
		   						        				supp_errors = 'Y',
		   						        				na.action = na.omit,
		   						        				lower = c(A = 0, tau_decay= 0))),
		   						       	tidied = map(fit, tidy),
					            		augment = map(fit, augment)

							)


					peekFit<-expFit%>%
			  			unnest(augment) 

		   			
			  		fitData<- peekFit
		        

			  		print('expDecayFit finished.')
	  	
		
		fitData	  
}

activeROIs<- n_transmission %>% dplyr::filter(n_total >= 2)#dplyr::filter(trackROI_key %in% saveROIs) 
ROI_list = unique(activeROIs$trackROI_key)

tmp_tag_var = c("G","H")	
      custom_y <- tribble(~sensor, ~normTime, ~dFF,
											  "GluSnFR3", 0, -0.1,
											  "GluSnFR3", 0, 1.2,
											  "JF646", 0, -0.1,
											  "JF646", 0.1, 0.7
                  								)
      						 custom_y$sensor = as.factor(as.character(custom_y$sensor))
						 levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))






print("Trying exp Decay Fit on JF646 peaks")
groupers = c('sensor','chemical_condition','trackROI_key','peak_key')
fittable_peaks <- traces_binned %>% dplyr::filter(trackROI_key %in% ROI_list)
print(fittable_peaks)
print(unique(fittable_peaks$sensor))
trace_binned_fitted <- expDecayFit(fittable_peaks)#, groupers)
trace_binned_fitted <- trace_binned_fitted[,!(names(trace_binned_fitted) %in% drops)]




subset_ROI_list<-  saveROIs


for(i in 1:length(subset_ROI_list)){
	current_ROI = subset_ROI_list[i]
	current_ROINumber = str_extract(current_ROI,"ROI\\d\\d\\d\\d")
	tmp_traces_binned = traces_binned %>% dplyr::filter(trackROI_key == current_ROI, chemical_condition == "Pre-Treatment")
	tmp_trace_binned_fitted = trace_binned_fitted %>% dplyr::filter(trackROI_key == current_ROI, chemical_condition == "Pre-Treatment")
	print(paste0("current ROI is : ", current_ROINumber))
	print(paste0("current chemical conditions are:", unique(tmp_trace_binned_fitted$chemical_condition)))
	print(paste0("current sensors available are:", unique(tmp_trace_binned_fitted$sensor)))
	# tmp_tau_decay_ms <- tmp_trace_binned_fitted %>% unnest(tidied) %>% 
	# 												dplyr::filter(term == "tau_decay") %>% 
	# 												select(sensor, chemical_condition, trackROI_key, peak_key, term, estimate) %>%
	# 												mutate(tau_decay_ms = estimate*1000) 
    # tmp_tau_decay_ms <- unique(tmp_tau_decay_ms) 
    # create_violin_plot(tmp_tau_decay_ms, "tau_decay_ms", expression(tau[decay]~"(ms)"), "#FF33FF", "X", paste0("violin__tau_decay_ms_",current_ROI))


		    avgtrace_Plot<-ggplot(tmp_traces_binned, aes(x=normTime, y = dFF,  group=peak_key, colour=sensor))+
					                                                    geom_line(alpha=0.3,size=1)+

					                                                    stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), geom="line", fun=mean, binwidth=interFrame.var*1.25,size=2,alpha=1)+
					                                                    #geom_path(data = tmp_trace_binned_fitted, aes(x=normTime, y = .fitted, group=peak_key),colour='black',linetype="dashed", size=2,alpha=0.2)+

					                                                    geom_blank(data=custom_y,aes(normTime,dFF), inherit.aes=FALSE)+															            
					                                                    
					                                                    labs( x="Normalized time (s)",
					                                                          y=expression(Delta*"F/F"),
					                                                          tag = tmp_tag_var[i],
					                                                          subtitle = paste0("Spine ",as.numeric(gsub("ROI00", "  ", current_ROINumber)))
					                                                          )+
					                                                          #title = unique(tmp_traces_binned$trackROI_key))+#,
					                                                          #subtitle = "Averaged Transmission Events")+
					                                                    scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-HTL-AM")))+
															            scale_x_continuous(expand=c(0,0),breaks=c(0,0.5),oob=scales::oob_keep)+ #limits=c(-0.4,0.6)
															            scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1))+#,limits=c(-0.1,0.6))+
															            coord_cartesian(clip='off',xlim=c(-0.2,0.7))+
															            
					                                                    theme_tufte()+
					                                                    my.theme+
					                                                    facet_grid(sensor~chemical_condition, labeller = label_parsed,scale="free_y")+
					                                                    theme(strip.text.y=element_blank(),
					                                                    	  legend.position = "none",
					                                                    	  strip.text.x = element_blank())

 					                                                    assign(paste0("avg_traces_dualcolor_",i), avgtrace_Plot, envir = .GlobalEnv)




 						save_plot(filename=paste0("averagePeaks_dualColor_clipped_",current_ROI),plot=avgtrace_Plot, plot_width=plot_dim*1.2,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)
					                         



 }





	get_tau_decay_ms <- trace_binned_fitted %>% unnest(tidied) %>% 
													dplyr::filter(term == "tau_decay") %>% 
													select(sensor, chemical_condition, trackROI_key, peak_key, term, estimate, dFF) %>%
													mutate(tau_decay_ms = estimate*1000,
															peak_color = "blue",
															max_dFF = max(dFF,na.rm=TRUE))
    get_tau_decay_ms <- get_tau_decay_ms %>% select(sensor, chemical_condition, trackROI_key, peak_key, term, estimate,tau_decay_ms, max_dFF)
    get_tau_decay_ms <- unique(get_tau_decay_ms) 
    my_medians_tau <- get_tau_decay_ms %>% group_by(chemical_condition) %>% summarise(median_tau_decay_ms = median(tau_decay_ms,na.rm=TRUE),
    																				  peak_color = "blue")



print("attempting to merge amplitude_dist_full with trace_binned_fitted")


# fix_amp_dist<-  amplitude_dist_full %>% dplyr::filter(chemical_condition == "Pre-Treatment", peak_color == "blue") %>% 
# 										group_by(chemical_condition, trackROI_key) %>% 
# 										mutate(index = row_number(interaction(chemical_condition,trackROI_key)),
# 												peak_key = paste("ctrl",trackROI_key, index, sep="-"),
# 												event_dt_ms = (time_of_JFPeak - time_of_GluPeak)*1000) %>%
# 										select(-index)	


# merge_amplitude_tau <- left_join(get_tau_decay_ms, fix_amp_dist) %>% dplyr::filter(chemical_condition == "Pre-Treatment")

tau_avgs = get_tau_decay_ms %>% dplyr::filter(chemical_condition == "Pre-Treatment") %>% group_by(chemical_condition, trackROI_key) %>% 
								summarise(n = n(),
											median_tau = median(tau_decay_ms, na.rm=TRUE),
											mean_tau = mean(tau_decay_ms, na.rm=TRUE),
											sd_tau = sd(tau_decay_ms, na.rm=TRUE),
											se_tau = sd_tau/sqrt(n))


subset_tau_decay<- get_tau_decay_ms %>% dplyr::filter(chemical_condition == "Pre-Treatment", trackROI_key %in% subset_ROI_list)


tau_labels<- c('Spine 4', 'Spine 20')
 tau_violin_plot<- ggplot(subset_tau_decay, aes(x = trackROI_key, y = tau_decay_ms)) +
							    geom_point(colour = "#FF33FF",size=4,alpha=1)+
							    stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
							    stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=6,fill="white", alpha=1, colour="black")+
                            	scale_x_discrete(labels = tau_labels)+
							    #geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
							    

							    labs(
							      y = expression(tau[decay]~"(ms)"),
							      x = "",
							      subtitle = "",#expression(JF[646]~tau[decay]),
							      tag = "I"
							    ) +
							    theme_tufte() +
							    my.theme +
							    theme(legend.position = "none",
							    	  axis.text.x=element_text(colour="black", size=scalefactor*34, family="sans") )

							  save_plot(filename = "tau_violin_plot_spine420_nice", plot = tau_violin_plot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

							  tau_violin <<- tau_violin_plot






	print(paste0("creating a histogram for: tau_decay_ms") )
								  upper_bound_x = 1500
								  lower_bound_x = 125
								  bin = (upper_bound_x - lower_bound_x) / 250
								  break_vec = seq(lower_bound_x,upper_bound_x, length.out = 5)
								  plot <- ggplot(get_tau_decay_ms, aes(x = tau_decay_ms, colour=chemical_condition)) +
								    geom_vline(data = my_medians_tau, aes(xintercept = median_tau_decay_ms, colour = chemical_condition), size = 2, lty = "solid") +
								    geom_freqpoly(aes(y = ..ncount.., colour=chemical_condition), binwidth = bin, alpha = 0.7, boundary = 0, size = 1) +
								    geom_histogram(aes(y = ..ncount.., fill=chemical_condition), binwidth = bin, alpha = 0.5, boundary = 0, position = "identity", colour = "black", size = 0.8) +
								    
								    labs(
								      x = expression(tau[decay]~"(ms)"),
								      y = expression(N/N[max]), #expression("Density"),#
								      fill = "",
								      tag = "X"#tag_var
								    ) +
								    scale_colour_manual(name = "", breaks=c("Pre-Treatment"),
				                                           		   values = c("Pre-Treatment"="#FF33FF"),
				                                           		   labels=c("Pre-Treatment" = "Fittable Peaks"),
				                                           		   guide="none")+
								    
								    scale_fill_manual(name = "", breaks=c("Pre-Treatment"),
				                                           		   values = c("Pre-Treatment"="#FF33FF"),
				                                           		   labels=c("Pre-Treatment" = "Fittable Peaks"),
				                                           		   guide="none")+
								    coord_cartesian(ylim = c(0, 1.3), xlim = c(lower_bound_x, upper_bound_x)) +
								    scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
								    scale_x_continuous(breaks=break_vec,limits =c(lower_bound_x,upper_bound_x))+
								    theme_tufte() +
								    #facet_grid(chemical_condition ~ ., labeller = label_parsed) +
								    my.theme +
								    theme(
								      legend.position = "top",
								      legend.direction = "vertical",
								      legend.justification = "center",
								      legend.title = element_blank(),
								      #axis.text.y = axis_text_y_arg,
								      #axis.title.y = axis_title_y_arg,
								      strip.text.y.right = element_blank()#element_text(colour="black", size=scalefactor*34, family="sans", angle=0)
								    )

								  save_plot(
								    filename = paste0("histo___tau_decay_ms", ""),
								    plot = plot,
								    plot_width=plot_dim*1.5,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)




    print(paste0("creating some additional avg trace plots based on different tau_decay_ms thresholds"))


    threshold_1 = 125
    threshold_2 = 2*threshold_1 #250
    threshold_3 = 3*threshold_1 #375
    threshold_4 = 4*threshold_1 #500
    threshold_5 = 5*threshold_1 #625
    threshold_6 = 7*threshold_1 #875

    tau_decay_1 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms <= threshold_1)
    tau_decay_2 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_1, tau_decay_ms <= threshold_2)
    tau_decay_3 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_2, tau_decay_ms <= threshold_3)
    tau_decay_4 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_3, tau_decay_ms <= threshold_4)
    tau_decay_5 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_4, tau_decay_ms <= threshold_5)
    tau_decay_6 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_5, tau_decay_ms <= threshold_6)
    tau_decay_7 <- get_tau_decay_ms %>% dplyr::filter(tau_decay_ms >= threshold_6)

    tracePlot_by_tau_decay<-function(traces_to_plot, tau_decay_thresh, filename,subtitle_expr){

    						peak_list<- unique(tau_decay_thresh$peak_key)
    						tmp_trace_df <- traces_to_plot %>% dplyr::filter(peak_key %in% peak_list)

    
      						
		    avgtrace_Plot<-ggplot(tmp_trace_df, aes(x=normTime, y = dFF,  group=peak_key, colour=sensor))+
					                                                    #geom_line(alpha=0.3,size=1)+

					                                                    stat_summary_bin(aes(group=interaction(sensor,chemical_condition)), geom="line", fun=mean, binwidth=interFrame.var*1.25,size=2,alpha=1)+
					                                                    
					                                                    labs( x="Normalized time (s)",
					                                                          y=expression(Delta*"F/F"),
					                                                          subtitle = subtitle_expr
					                                                          #tag = tmp_tag_var[i],
					                                                          #subtitle = paste0("Spine ",as.numeric(gsub("ROI00", "  ", current_ROINumber)))
					                                                          )+
					                                                          #title = unique(tmp_traces_binned$trackROI_key))+#,
					                                                          #subtitle = "Averaged Transmission Events")+
					                                                    scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-HTL-AM")))+
															            scale_x_continuous(expand=c(0,0),breaks=c(0,0.5),oob=scales::oob_keep)+ #limits=c(-0.4,0.6)
															            scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1))+#,limits=c(-0.1,0.6))+
															            coord_cartesian(clip='off',xlim=c(-0.2,0.7),ylim=c(-0.2,0.7))+
															            
					                                                    theme_tufte()+
					                                                    my.theme+
					                                                    facet_grid(sensor~chemical_condition, labeller = label_parsed,scale="free_y")+
					                                                    theme(strip.text.y=element_blank(),
					                                                    	  legend.position = "none",
					                                                    	  strip.text.x = element_blank())

							save_plot(filename=paste0(filename),plot=avgtrace_Plot, plot_width=plot_dim,plot_height=plot_dim*2, scale_f = scalefactor,dpi_val = 300)


     }

tracePlot_by_tau_decay(traces_binned, tau_decay_1, "tau_decay-125ms-avg_traces",subtitle_expr = expression(tau[decay]~"< 125 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_2, "tau_decay-125ms_250ms-avg_traces",subtitle_expr = expression(tau[decay]~"125-250 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_3, "tau_decay-250ms_375ms-avg_traces",subtitle_expr = expression(tau[decay]~"250-375 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_4, "tau_decay-375ms_500ms-avg_traces",subtitle_expr = expression(tau[decay]~"375-500 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_5, "tau_decay-500ms_625ms-avg_traces",subtitle_expr = expression(tau[decay]~"500-625 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_6, "tau_decay-625ms_875ms-avg_traces",subtitle_expr = expression(tau[decay]~"625-875 ms"))
tracePlot_by_tau_decay(traces_binned, tau_decay_7, "tau_decay-875ms+-avg_traces",subtitle_expr = expression(tau[decay]~"> 875 ms"))





## pseudo-code to generate extended data Figure 6

# ROI_avgs should include corr-slope for both 
# instantiate a separate for loop for plotting the before and after treatment.
# loop 1: vehicle
	png_dir = "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_vFinal/synTransmission_v1/forMapping/_registered_11"
	tmp_tag_var = c("A", "B")
	var_to_plot = c('corr_slope', 'corr_slope')
	
# loop 2: APV
	png_dir = "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_vFinal/synTransmission_v1/forMapping/_registered_6"
	tmp_tag_var = c("C", "D")
	var_to_plot = c('corr_slope', 'corr_slope')



# c("A","B","C","D") 


#
		ext_data_saveROIs<- c("dish8-plate07-region05-ROI0020",
									"dish7-plate10-region04-ROI0016")


					#"Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v6/synTransmission_v2/forMapping/_registered_3"
		png_dirs = c("Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v6/synTransmission_v2/forMapping/_registered_10",
            		 "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v6/synTransmission_v2/forMapping/_registered_6")


        prefix_str = "AVG_GluSnFR3_JF646_4Ca_cis_"
        suffix_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        all_lut_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        dir_regex = "_registered_(.*?)"
        save_bool=TRUE
         
        
        
        

        ###scale bar params! set the left-corner 
        user_x = 51.2-10.6
        user_y = 50.2


        #tmp_tag_var = c("A","B","C","D")
        #vars_to_plot<- c('corr_slope', "corr_slope","corr_slope","corr_slope")#c("releaseProbability_perROI")
        #guide_limits = list(c(0,1500),c(0,1500),c(0,1500),c(0,1500))#,#c(0,0.25),c(0,1))
        #guide_titles = c(expression("Signal Correlation"),expression("Signal Correlation"),expression("Signal Correlation"),expression("Signal Correlation"))



           print("Taking a crack at generating png plots for extended data Figure 6")
     

        for( j in 1:length(png_dirs)){
        	current_png_dir = png_dirs[j]
        	print(paste0("Current png_dir is : ", current_png_dir))

        	if(j == 1){
        		tmp_tag_var = c("A","B")
        		vars_to_plot<- c('corr_slope', "corr_slope")
        		guide_limits = list(c(0,1.0),c(0,1.0))
        		guide_titles = c(expression(atop("Correlation","slope")),expression(atop("Correlation","slope")))
        		ROI_to_be_saved = ext_data_saveROIs[j]
        		print(paste0("ROI to be saved is: ",ROI_to_be_saved))
        		ROInames <- str_extract(ROI_to_be_saved,"ROI\\d\\d\\d\\d")
        		vid_keys_to_save <- unique(str_extract(ROI_to_be_saved, "dish(.*?)region\\d\\d") )
        


        		subset_df<- ROI_avgs %>% mutate(vid_key =str_extract(trackROI_key, "dish(.*?)region\\d\\d"),
        										ROINumber = str_extract(trackROI_key, "ROI\\d\\d\\d\\d") ) %>% 
        								 dplyr::filter(vid_key %in% vid_keys_to_save)
        		vid_key_vec <- unique(subset_df$vid_key)
        		print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))

        

        	}


        	if(j == 2){
        		tmp_tag_var = c("E","F")
        		vars_to_plot<- c('corr_slope', "corr_slope")
        		guide_limits = list(c(0,1.0),c(0,1.0))
        		guide_titles = c(expression(atop("Correlation","slope")),expression(atop("Correlation","slope")))
        		ROI_to_be_saved = ext_data_saveROIs[j]
        		ROInames <- str_extract(ROI_to_be_saved,"ROI\\d\\d\\d\\d")
        		vid_keys_to_save <- unique(str_extract(ROI_to_be_saved, "dish(.*?)region\\d\\d") )
				

				subset_df<- ROI_avgs %>% mutate(vid_key =str_extract(trackROI_key, "dish(.*?)region\\d\\d"),
        										ROINumber = str_extract(trackROI_key, "ROI\\d\\d\\d\\d") ) %>% 
        								 dplyr::filter(vid_key %in% vid_keys_to_save)
				vid_key_vec <- unique(subset_df$vid_key)
                print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))



	        	}
        

	        for(k in 1:length(vars_to_plot)){      


                                #guide_bool = TRUE
                                var_to_plot = vars_to_plot[k]
                                tmp_guide_title = guide_titles[k]
                                tmp_limits = guide_limits[[k]]
                                tmp_tag = tmp_tag_var[k]
                                legend.pos = c(0.13, 0.55)

                                if(k == 1){
                                	tmp_subset_transmission<- subset_df %>% dplyr::filter(chemical_condition == "Pre-Treatment")
                                	chem_key = "pre_treat"
                                }
                                if(k == 2){
                                	tmp_subset_transmission<- subset_df %>% dplyr::filter(chemical_condition != "Pre-Treatment")
                                	chem_key = "wash-in"
                                }
                                if(k == 1){tmp_scalebar = TRUE
                                			guide_bool = TRUE
                                			} else {tmp_scalebar = FALSE
                                					guide_bool = FALSE
                                				     }
                                tmp_plot = png_plotter(df=tmp_subset_transmission,
                                                                    png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                                                    vars_to_plot=var_to_plot, 
                                                                    save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
                                                                    ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar,add_arrows=TRUE,
                                                                    scalebar_length=10,binfactor=2,chem_key = chem_key, legend.pos = legend.pos)



                        if(j == 1 & k == 1){
                            region_vehicle_A  <<- tmp_plot
                            }
                        if(j == 1 & k == 2){
                            region_vehicle_B  <<- tmp_plot
                            }

                   		if(j == 2 & k == 1){
                            region_APV_E  <<- tmp_plot
                            }
                        if(j == 2 & k == 2){
                            region_APV_F  <<- tmp_plot
                            }



             }
        }


print("finished attempting Extended Figure 6")









##### GET SPATIAL DATA ##### ROI OVERLAY + ARROW

# ###### for loop for plotting spatial data
        png_dir = "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v6/synTransmission_v2/forMapping/_registered_3"#,
                        #"Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_6",
                        #"Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_11")
        prefix_str = "AVG_GluSnFR3_JF646_4Ca_cis_"
        suffix_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        all_lut_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        dir_regex = "_registered_(.*?)"
        save_bool=TRUE
         
        print("Taking a crack at generating png plots")
        print(saveROIs)
        vid_keys_to_save <- unique(str_extract(saveROIs, "dish(.*?)region\\d\\d") )
        print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))
        
        ROI_avg_clean <- ROI_avgs %>% dplyr::filter(chemical_condition == "Pre-Treatment")
        tau_avg_clean <- tau_avgs %>% ungroup() %>% select(trackROI_key,median_tau,mean_tau,se_tau)
        ROI_avg_merge <- left_join( ROI_avg_clean, tau_avg_clean)

        subset_df<- ROI_avg_merge %>% mutate(vid_key =str_extract(trackROI_key, "dish(.*?)region\\d\\d"),
        								ROINumber = str_extract(trackROI_key, "ROI\\d\\d\\d\\d") ) %>% 
        						dplyr::filter(vid_key %in% vid_keys_to_save,chemical_condition == "Pre-Treatment")
        ROInames <- str_extract(saveROIs,"ROI\\d\\d\\d\\d")
        print(ROInames)
        vid_key_vec <- unique(subset_df$vid_key)
        vars_to_plot<- c('mean_Glu_dFF', "mean_JF_dFF","median_tau","avg_n","corr_slope")#c("releaseProbability_perROI")"corr_slope",
        save_bool=TRUE
        guide_bool = TRUE

        guide_limits = list(c(0,1),c(0,0.5),c(0,500),c(0,10),c(0,1.5))#,#c(0,0.25),c(0,1))c(0,1.5),
        guide_titles = c(expression(iGlu~Delta*"F/F"),expression(JF[646]~Delta*"F/F"),expression(tau[decay]*", ms"),expression(N[transmission]),expression("Signal Correlation"))#,,expression(JF[646]~(Hz)),expression(P[transmission]))# # expression(CV[tau*"decay"])#
        
        #expression("Avg."~Evoked[Delta*"F/F"] / Spont[Delta*"F/F"])
        tmp_tag_var = c("","F","F","F","")


        ###scale bar params! set the left-corner 
        user_x = 51.2-10.6
        user_y = 50.2
        #tag_vars = c("F","G"    )



        #tag_vars = c("H","F")
        for(j in 1:length(png_dir)){
            current_png_dir = png_dir
            print(paste0("current png directory is: ",current_png_dir) )
            print(paste0("Vid_key is : ",vid_key_vec[j]) )
               
            tmp_subset_transmission <- subset_df %>% dplyr::filter(vid_key == vid_key_vec[j]) %>% mutate(fix_Ca = "4Ca")

            print(head(tmp_subset_transmission))

            for(k in 1:length(vars_to_plot)){      


                                guide_bool = TRUE
                                var_to_plot = vars_to_plot[k]
                                tmp_guide_title = guide_titles[k]
                                tmp_limits = guide_limits[[k]]
                                tmp_tag = tmp_tag_var[k]
                                print(paste0("tmp_tag_var = ", tmp_tag))
                                if(k == 1){tmp_scalebar = TRUE} else {tmp_scalebar = FALSE}
                                tmp_plot = png_plotter(df=tmp_subset_transmission,
                                                                    png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                                                    vars_to_plot=var_to_plot, 
                                                                    save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
                                                                    ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar,add_arrows=TRUE,
                                                                    scalebar_length=10,binfactor=2)



                        if(k == 4){
#### save to global env
                            regionF  <<- tmp_plot


                            }
#                         if(k == 2){
# #### save to global env
#                            regionG  <<- tmp_plot

#                             }
                        # if(k == 3){
                        # 	regionH <<- tmp_plot
                        # }         



             }
        }


						


						#output_list = list(avg_transmission_stats = my_avgs_full, ROI_avgs = ROI_avgs, ROI_avgs_pretreat = ROI_avgs_pretreat)
						#output_list
						#traces_binned
                        #output_list = list(peak_key_maxFinder = amplitude_dist_full)
                        #output_list
                        output_list = list(avg_transmission_stats = my_avgs_full, 
                        					ROI_freq = frequency_dist, 
                        					trace_times = trace_times, 
                        					amplitude_dist_full = amplitude_dist_full, 
                        					n_transmission = n_transmission,
                        					ROI_avgs = ROI_avgs, 
                        					traces_binned = traces_binned, 
                        					traces_binned_fitted = trace_binned_fitted, 
                        					get_tau_decay_ms = get_tau_decay_ms, 
                        					#merge_amplitude_tau = merge_amplitude_tau,
                        					#fix_amp_dist = fix_amp_dist,
                        					#tau_avgs = tau_avgs,
                        					ROI_avg_merge = ROI_avg_merge)
                        output_list

					}



transmission_output<-transmission_stats(df_traces, saveROIs = ROIs_to_save)





