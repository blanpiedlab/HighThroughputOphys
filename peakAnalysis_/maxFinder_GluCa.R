#### maxFinder_GluCa_v1.R 

# Written 5.14.2024 by Samuel T Barlow
# 
# The purpose of this script is to define a function which identifies maxima in JF646 recordings based on the position (in time) of iGluSnFR3 signals identified by peakFinder.R. 
# This function builds on single_ROI_display_v2.R, which uses a similar approach to perform single dendritic spine analysis.
# However, in maxFinder_GluCa_v1.R, we want to build on this approach to make it scalable over the entire dataset. 

# maxFinder_GluCa_v1 should take two inputs: 
# 1) df, a dataframe containing all correlated 2-channel recordings. 
# 2) window_boundary, a vector ( c(start, end)) that defines the temporal boundaries around iGluSnFR3 peaks in which to search for a JF646 maximum. Use time in seconds. Base implementation is c(0.01, 0.25), referring to 10 ms prior and 250 ms after a iGluSnFR3 peak.

# maxFinder_GluCa_v1 should have multiple outputs: 
# 1) df_maxima, a dataframe featuring all of the found maxima positions.
# 2) df_traces_annotated, a dataframe featuring all of the traces + peak_color, a few other columns
# 3) window_info, a dataframe featuring all of the calculated window information, which should be output for quality control purposes. (Ca_plot_elements) 

# to verify implementation, there will be a toggle-able output (Boolean) which uses single_ROI_display_v2 code to graph the output. This will be identical to the output for single_ROI_display_v2.R, so this should only need to be run during QC. 


# approach: 
# maxFinder_GluCa_v1 will use the existing approach from single_ROI_display_v2.R, which relies on a for loop implementation to break down the dataset into the following iterable parts: 
# ROIs - Each trackROI_key will be analyzed one at a time to ensure that there is no cross-talk between ROIs when defining peak maxima. 
# chemical_condition - Most trackROI_keys in the dataset feature a "Pre-treatment" condition before manipulation, and a post-treatment condition (vehicle or APV). To avoid cross-talk between data streams, each trackROI_key will have individual chemical conditions analyzed separately. 


source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))



library(ggforce)
library(scales)
library(ggthemes)
library(ggpmisc)
library(ggpubr)




maxFinder<- function(df, window_boundary = c(0.01,0.25), graph_output = FALSE){

			sd_threshold = 3.5			
			ROI_list<- unique(df$trackROI_key)

			total_ROIs<- length(ROI_list)

			df_maxima<- list()
			df_traces_annotated<- list()
			window_info<- list()

			for(i in 1:total_ROIs){
									current_ROI <- ROI_list[i]


									tmp_trace = df %>% dplyr::filter(trackROI_key %in% current_ROI)

									#report which ROI we are working on
									graphTitle = first(unique(tmp_trace$trackROI_key) )
									print(paste0("Running maxFinder_GluCa for ROI : ", graphTitle))


									#special case for broken video to avoid off-target info
									if(unique(tmp_trace$vid_key) == 'dish7-plate07-region04'){

							                            	tmp_trace_brokenROI<- tmp_trace %>% dplyr::filter(chemical_condition == "APV", vid_key == "dish7-plate07-region04", absoluteTime >= 6) ###this dataset had a focus issue that is messing with quantification
											                tmp_trace_ctrl <- tmp_trace %>% dplyr::filter(chemical_condition == "ctrl",vid_key == "dish7-plate07-region04")
											                tmp_trace_fixed <- bind_rows(tmp_trace_brokenROI,tmp_trace_ctrl)
			                            					tmp_trace <- tmp_trace_fixed
			                            					rm(tmp_trace_brokenROI,tmp_trace_ctrl,tmp_trace_fixed)
							                            }

							         get_baseline_indices <- tmp_trace %>% dplyr::filter(sensor == "JF646",alpha_str == "notPeak")
									 baseline_times <- get_baseline_indices %>% ungroup() %>% select(chemical_condition, trackROI_key, absoluteTime)
									 baseline_pairs<-suppressMessages(inner_join(tmp_trace,baseline_times) )
									 get_hline<- baseline_pairs %>% dplyr::filter(sensor == "JF646", chemical_condition == "ctrl") %>% summarise(hline_mean = mean(dFF,na.rm=TRUE),
																																				 hline_sd = sd_threshold*sd(dFF,na.rm=TRUE),
		 															        																	 hline_tresh = hline_mean+hline_sd)
									 hline_val <- get_hline$hline_tresh
									 hline_mean<- get_hline$hline_mean			



									 ##### This is the part of the code that is actually finding maxima in a defined window. 

									 def_Glu_window_boundaries <- tmp_trace %>% 
																			dplyr::filter(peakID != "NotPeak") %>% # only peaks as defined by peakFinder.R
			 															    dplyr::filter(sensor == "GluSnFR3") %>% # only iGluSnFR3 traces to define the temporal boundaries of the trace in which to look for JF646 maxima.
																	        group_by(sensor,trackROI_key,chemical_condition,peakID) %>% # individualize the operation to single peakIDs
																	        slice(which.max(dFF)) %>% #pull out just the time of peak maximum
																	        ungroup() %>%  
																	        group_by(chemical_condition,trackROI_key) %>% 
																	        arrange(absoluteTime) %>% # probably unnecessary, but orders the dataframe by absoluteTime within ROIs.  
																	        mutate(check_time_start = absoluteTime-window_boundary[1],
																	        	   check_time_end = absoluteTime + window_boundary[2]) %>%
																	        ungroup() %>%
																	        select(sensor, chemical_condition,trackROI_key,absoluteTime,dFF,check_time_start,check_time_end) 
									 
				                     # apply iGluSnFR3 times to Ca_trace
									 Ca_trace <- tmp_trace %>% dplyr::filter(sensor == "JF646")

									 start_times  = def_Glu_window_boundaries %>% 
									 										ungroup() %>% 
																	        group_by(chemical_condition, trackROI_key)  %>% 
																	        select(chemical_condition, trackROI_key,check_time_start) %>%
																	        distinct(check_time_start)
		                             end_times = def_Glu_window_boundaries %>% 
									                                        ungroup() %>% 
																	        group_by(chemical_condition, trackROI_key)  %>% 
																	        select(chemical_condition, trackROI_key,check_time_end) %>%
																	        distinct(check_time_end)
						             n_intervals = nrow(start_times)
									 #instantiate a list for the interval data
									 selected_intervals <- list()
									 

									 ## use for loop (bleh) to define bins of data.frame
									 for (j in 1:n_intervals) {
									 							current_chem_condition = start_times$chemical_condition[j]
																bin_number = j
																bin_start = start_times$check_time_start[j]
																bin_end = end_times$check_time_end[j]
																
																selected_data <- Ca_trace %>% ungroup() %>% 
																							  dplyr::filter(chemical_condition == current_chem_condition, absoluteTime >= bin_start & absoluteTime <= bin_end) %>%
																	                    	  mutate(window_start = first(absoluteTime),
																	                    	         window_end = last(absoluteTime),
																	                    			 bin = bin_number)
																	                    			 #peak_key = paste(chemical_condition,trackROI_key,bin,sep="-"))

		 														selected_intervals[[j]] <- selected_data

			                         }

			                         Ca_trace_binned <- bind_rows(selected_intervals) 
			                         Ca_plot_elements<- Ca_trace_binned %>% ungroup() %>% select(sensor, chemical_condition, trackROI_key,window_start, window_end)
									 Ca_plot_elements<-unique(Ca_plot_elements)

									 get_JF_peakPositions <- Ca_trace_binned %>% 
																			group_by(sensor,chemical_condition,trackROI_key,bin) %>%
																	        slice(which.max(dFF)) %>%
																	        mutate(peak_color = ifelse(dFF >= hline_val, "blue","white"))




				                   	               	

				                	 #define baseline dataframe
				                   	 Glu_baseline<- baseline_pairs %>% dplyr::filter(sensor == "GluSnFR3") %>% ungroup() %>% select(chemical_condition, trackROI_key,absoluteTime,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
				                     colnames(Glu_baseline)[colnames(Glu_baseline) == "dFF"] <- "Glu_dFF"
				                     colnames(Glu_baseline)[colnames(Glu_baseline) == "absoluteTime"] <- "time_of_GluPeak"
				                     
				                     JF_baseline<- baseline_pairs %>% dplyr::filter(sensor == "JF646") %>% ungroup() %>% select(chemical_condition, trackROI_key,absoluteTime,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
				                     colnames(JF_baseline)[colnames(JF_baseline) == "dFF"] <- "JF_dFF"
				                     colnames(JF_baseline)[colnames(JF_baseline) == "absoluteTime"] <- "time_of_JFPeak"

				                     check_baseline_dist<- suppressMessages(left_join(Glu_baseline,JF_baseline) %>% mutate(peak_color = "baseline") )


				                     ## define peak position dataframe
				                     get_Glu_peakPositions = def_Glu_window_boundaries              


				                     Glu_peaks <- get_Glu_peakPositions %>% ungroup() %>% select(chemical_condition, trackROI_key,absoluteTime,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
				                     colnames(Glu_peaks)[colnames(Glu_peaks) == "dFF"] <- "Glu_dFF"
				                     colnames(Glu_peaks)[colnames(Glu_peaks) == "absoluteTime"] <- "time_of_GluPeak"
				                     JF_peaks <- get_JF_peakPositions %>% ungroup() %>% select(chemical_condition, trackROI_key,absoluteTime,dFF,peak_color) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
				                     colnames(JF_peaks)[colnames(JF_peaks) == "dFF"] <- "JF_dFF"
				                     colnames(JF_peaks)[colnames(JF_peaks) == "absoluteTime"] <- "time_of_JFPeak"

				                     amplitude_dist = suppressMessages(left_join(Glu_peaks,JF_peaks, by=c("index","chemical_condition","trackROI_key")) )
				                     amplitude_dist = bind_rows(amplitude_dist,check_baseline_dist)
				                     amplitude_dist$chemical_condition = factor(amplitude_dist$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
									 
				                     amplitude_dist_summary<- amplitude_dist %>% dplyr::filter(peak_color == "blue") %>% 
				                     											 group_by(trackROI_key) %>% 
				                     											 summarise(mean_amplitude_JF = mean(JF_dFF,na.rm=TRUE),
				                     											 			sd_amplitude_JF = sd(JF_dFF, na.rm=TRUE),
				                     											 			
				                     											 			mean_amplitude_Glu = mean(Glu_dFF,na.rm=TRUE),
				                     											 			sd_amplitude_Glu = sd(Glu_dFF, na.rm=TRUE)
				                     											 			)


									amplitude_dist_z_score<- left_join(amplitude_dist, amplitude_dist_summary) %>% group_by(trackROI_key, chemical_condition) %>% mutate(z_score_JF = (JF_dFF - mean_amplitude_JF)/sd_amplitude_JF,
																																												z_score_Glu = (Glu_dFF - mean_amplitude_Glu)/sd_amplitude_Glu)
																																										 		
			                     

									 ### This code block is copy-pasted from singleROI_display_v2.R ###
				                    if(graph_output == TRUE){

				                    		subDir <- "test_maxFinder_output_v1"   
											source(paste0(path,"sub-functions/setFigureDirectory.R")) 

				                    		scalefactor = 0.75
											plot_dim = 16

										 	strip_text_arg_y = element_blank()
										 	strip_text_arg_x = element_text(colour="black", size=34*scalefactor, family="sans")
										 	legend_arg = "none"
										 	break_arg = c(0,1,2,3)


											amplitude_scatterplot<-ggplot(amplitude_dist, aes(x=Glu_dFF, y=JF_dFF, group=interaction(trackROI_key,chemical_condition),colour=peak_color,alpha=peak_color))+
																		  	geom_point(size=4)+
																		  	geom_smooth(data=subset(amplitude_dist,peak_color == "blue"), method = "lm",formula = y~ x, se = FALSE, colour = "black", size=1.5,alpha=0.8)+
																		  	geom_hline(yintercept =hline_val,lty="dashed",size=0.5)+
																		    labs( 	x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
																		    	    y=expression(JF[646]*" amplitude ("*Delta*"F/F)"),
																		            tag="B")+
																		    scale_colour_manual(name = "", 
																		    					breaks=c("blue","white", "baseline"),
																		    					values = c("blue"='blue',"white" = "grey71","baseline" = "grey71"),
																		                        labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"))+
																		    scale_alpha_manual(name = "", 
																		    					values = c("blue"=0.9,"white" = 0.6,"baseline" = 0.4),
																		    					labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), 
																		    					guide="none")+																	                                           
																		    coord_cartesian(xlim=c(0,2.5),ylim=c(-0.1,1.0))+
																		    scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
																		    scale_x_continuous(breaks=break_arg)+
																		    facet_grid(~chemical_condition,labeller = label_parsed)+
																		    theme_tufte()+
																		    my.theme +
																		    theme(legend.position = legend_arg,
																		    		strip.text.x = strip_text_arg_x)








											z_score_scatterplot<-ggplot(amplitude_dist_z_score, aes(x=z_score_Glu, y=z_score_JF, group=interaction(trackROI_key,chemical_condition),colour=peak_color,alpha=peak_color))+
																		  	geom_point(size=4)+
																		  	geom_smooth(data=subset(amplitude_dist_z_score,peak_color == "blue"), method = "lm",formula = y~ x, se = FALSE, colour = "black", size=1.5,alpha=0.8)+
																		  	#geom_hline(yintercept =hline_val,lty="dashed",size=0.5)+
																		    labs( 	x=expression("iGluSnFR3 signal z-score"),
																		    	    y=expression(JF[646]*"signal z-score"),
																		            tag="C")+
																		    scale_colour_manual(name = "", 
																		    					breaks=c("blue","white", "baseline"),
																		    					values = c("blue"='blue',"white" = "grey71","baseline" = "grey71"),
																		                        labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"))+
																		    scale_alpha_manual(name = "", 
																		    					values = c("blue"=0.9,"white" = 0.8,"baseline" = 0.8),
																		    					labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), 
																		    					guide="none")+																	                                           
																		   # coord_cartesian(xlim=c(0,2.5),ylim=c(-0.1,1.0))+
																		   # scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
																		    #scale_x_continuous(breaks=break_arg)+
																		    facet_grid(~chemical_condition,labeller = label_parsed)+
																		    theme_tufte()+
																		    my.theme +
																		    theme(legend.position = legend_arg,
																		    		strip.text.x = strip_text_arg_x)
								            min_val = -0.1
																
											facet.labels <- c("iGluSnFR3","JF646-BAPTA-AM")
											names(facet.labels)<- c("GluSnFR3", "JF646")
											tmp_trace$chemical_condition = factor(tmp_trace$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
											Ca_plot_elements$chemical_condition = factor(Ca_plot_elements$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
											get_Glu_peakPositions$chemical_condition = factor(get_Glu_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
											get_JF_peakPositions$chemical_condition = factor(get_JF_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
											hline_df<- data.frame(sensor = "JF646", hline = hline_val, min = min_val)
											# get_times<- tmp_trace %>% dplyr::filter(sensor == 'JF646') %>% ungroup() %>% select(sensor,chemical_condition,trackROI_key,absoluteTime) %>% 
											# 																					group_by(sensor, chemical_condition,trackROI_key) %>%
											# 																					summarise(max_time = max(absoluteTime,na.rm=TRUE))
											# hline_df<- left_join(get_times,hline_df)





										 	custom_y <- tribble(~sensor, ~absoluteTime, ~dFF,
										  						"GluSnFR3", 0, min_val,
										  						"GluSnFR3", 0, 1.2,
										  						"JF646", 0, min_val,
										  						"JF646", 0, 0.75
										                    	)

											my_y_breaks <- function(x) {if(max(x) <= 1.1)  seq(0.0,1.0,by=0.2) else if(max(x) <= 2 & max(x) >= 1.2) seq(0.0,2,by=0.5) else seq(0.0,7.5,by=2.5)}
											my_breaks <- function(x) { if(max(x) < 20)  seq(0,15,5) else seq(0,25,5) }
											formatter <- function(...){function(x) format(round(x, 1), ...)  }

											tmp_trace$sensor = as.factor(as.character(tmp_trace$sensor))
											levels(tmp_trace$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))
																 
											Ca_plot_elements$sensor = as.factor(as.character(Ca_plot_elements$sensor))
											levels(Ca_plot_elements$sensor) <- c(expression(JF[646]*"-BAPTA-HTL-AM"))

											get_Glu_peakPositions$sensor = as.factor(as.character(get_Glu_peakPositions$sensor))
											levels(get_Glu_peakPositions$sensor) <- c("iGluSnFR3")

											get_JF_peakPositions$sensor = as.factor(as.character(get_JF_peakPositions$sensor))
											levels(get_JF_peakPositions$sensor) <- c(expression(JF[646]*"-BAPTA-HTL-AM"))


											custom_y$sensor = as.factor(as.character(custom_y$sensor))
											levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-HTL-AM"))

											hline_df$sensor = as.factor(as.character(hline_df$sensor))
											levels(hline_df$sensor) <- c(expression(JF[646]*"-BAPTA-HTL-AM"))
																 


											tracePlot<-ggplot(tmp_trace, aes(x=absoluteTime, y = dFF, group = sensor)) +
																            geom_rect(data=Ca_plot_elements, aes(xmin = window_start, xmax = window_end,group=chemical_condition), ymin = -Inf, ymax = Inf,alpha=0.14, fill="grey21",inherit.aes=FALSE)+
																            geom_vline(data=Ca_plot_elements,aes(xintercept = window_start,group=chemical_condition), colour = "black", alpha = 0.7, linetype = "longdash", size =0.2) +
																            geom_vline(data=Ca_plot_elements,aes(xintercept = window_end,group=chemical_condition), colour = "black", alpha = 0.7, linetype = "longdash", size =0.2) + 
																            geom_blank(data=custom_y,aes(absoluteTime,dFF))+
																            geom_line(aes(colour = sensor),size = 1.2) + 
																           
																            geom_point(data=get_Glu_peakPositions, aes(x=absoluteTime,y=dFF), shape=23, fill="blue",alpha=0.7,colour="black", size=5)+
																            geom_point(data=get_JF_peakPositions, aes(x=absoluteTime,y=dFF,fill=peak_color,shape=peak_color), colour="black",alpha=0.7, size=5)+ #colour="blue",fill="blue", shape=23,
																            
																            labs( x="Time (s)",
																                  y=expression(Delta*"F/F"), 
																                  tag = "A")+
																                  
																            scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-HTL-AM")))+
																            scale_fill_manual(values = c("blue","grey71"),labels=c("blue","white"))+
																            scale_shape_manual(values = c(23,23),labels=c("blue","white"))+
																            scale_x_continuous(expand=c(0,0))+
																            coord_cartesian(clip='off')+
																            theme_tufte()+
																            my.theme+
																            
																            theme(legend.position=legend_arg,
																             		strip.text.x=strip_text_arg_x,
																             		strip.text.y= strip_text_arg_y)+
																            facet_grid(sensor~chemical_condition, labeller = label_parsed,scales="free_y")+
																            scale_y_continuous(labels = formatter(nsmall=1), breaks= my_y_breaks) 


																


											gs = list(tracePlot,  #1
						                                 amplitude_scatterplot,
						                                 	z_score_scatterplot)  #3
						                    margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))              
						                    hlay <- rbind(c(1,1,1,1,1,1,1,2,2,2,2),
						                                       c(1,1,1,1,1,1,1,2,2,2,2),
						                                       c(1,1,1,1,1,1,1,2,2,2,2),
						                                       c(1,1,1,1,1,1,1,2,2,2,2),
						                                       c(3,3,3,3,3,NA,NA,NA,NA,NA,NA),
						                                       c(3,3,3,3,3,NA,NA,NA,NA,NA,NA),
						                                       c(3,3,3,3,3,NA,NA,NA,NA,NA,NA),
						                                       c(3,3,3,3,3,NA,NA,NA,NA,NA,NA)
						                                        

						                                   )
						                
						                    

						                    
						                    ROI_report = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
						                    save_plot(filename=paste0("ROI_report_",graphTitle,"_"),plot=ROI_report, plot_width=plot_dim*3,plot_height=plot_dim*2, scale_f = scalefactor,dpi_val = 150)
						                    print(paste0("Should be done with saving out plots for : ", graphTitle ) )

											}

									 #At this point, we now have a dataframe of the Glu_maxima vs. JF_maxima for individual peaks and baseline points for a single ROI for each chemical condition tested. (df_maxima[i])
			                     # amplitude_dist should be output into a list for each iteration of the script. 

								 df_maxima[[i]] <- amplitude_dist_z_score
								 #df_traces_annotated[[i]]<- tmp_trace #This is a bit extraneous computationally, since the underlying trace doesn't get modified in this version of the script.   
								 window_info[[i]] <- Ca_plot_elements

								 #rm(list=setdiff(ls(), c("df_maxima",'window_info',"ROI_list","total_ROIs")))  #clear all vars in memory except for flagged data.

			}


		                     
			df_maxima_final <- bind_rows(df_maxima)
			window_info_final<- bind_rows(window_info)

			output_data<- list(df_maxima_final, window_info_final)
			output_data         

}
