#### peak_clipper.R 

# Written 5.14.2024 by Samuel T Barlow
# 
# The purpose of this script is to build a dataframe that has clipped out a window around iGluSnFR3 peaks and visualized the same region in corresponding JF646 traces. We will use this clipping strategy to normalize time for every peak in the dataset, align these and generate an averaged trace.   


# peak_clipper should take two inputs: 
# 1) df, a dataframe containing all correlated 2-channel recordings. 
# 2) window_boundary, a vector ( c(start, end)) that defines the temporal boundaries around iGluSnFR3 peaks in which to search for a JF646 maximum. Use time in seconds. Base implementation is c(0.100, 0.5), referring to 100 ms prior and 500 ms after a iGluSnFR3 peak.

# peak_clipper should have one output : 
# 1) df_traces_annotated, a dataframe featuring only the clipped transients now normalized in time 
# to verify implementation, there will be a toggle-able output (Boolean) which uses single_ROI_display_v2 code to graph the output. 




peak_clipper<- function(df, window_boundary = c(0.1,0.5)){

			
			ROI_list<- unique(df$trackROI_key)

			total_ROIs<- length(ROI_list)

			#df_maxima<- list()
			df_traces_annotated<- list()
			#window_info<- list()

			for(i in 1:total_ROIs){
									current_ROI <- ROI_list[i]


									tmp_trace = df %>% dplyr::filter(trackROI_key %in% current_ROI)

									#report which ROI we are working on
									graphTitle = first(unique(tmp_trace$trackROI_key) )
									print(paste0("Loop #",i,": running peak_clipper for ", graphTitle))


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
																																				 hline_sd = 3.5*sd(dFF,na.rm=TRUE),
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
									 Glu_trace <- tmp_trace %>% dplyr::filter(sensor == "GluSnFR3")
									 Ca_trace <- tmp_trace %>% dplyr::filter(sensor == "JF646") %>% ungroup() %>% group_by(chemical_condition,trackROI_key)
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
									 #selected_Ca_intervals <- list()
									 

									 ## use for loop (bleh) to define bins of data.frame
									 for (j in 1:n_intervals) {
									 							current_chem_condition = start_times$chemical_condition[j]
																bin_number = j
																bin_start = start_times$check_time_start[j]
																bin_end = end_times$check_time_end[j]
																
																selected_Glu_data <- Glu_trace %>% ungroup() %>%
																							  dplyr::filter(chemical_condition == current_chem_condition, absoluteTime >= bin_start & absoluteTime <= bin_end) %>%
																							  group_by(chemical_condition,trackROI_key) %>%
																	                    	  mutate(window_start = first(absoluteTime),
																	                    	         window_end = last(absoluteTime),
																	                    			 bin = bin_number,
																	                    			 time_of_max = absoluteTime[which.max(dFF)],
																	                    			 max_dFF = max(dFF,na.rm=TRUE),
																	                    			 #peak_color = "blue",
																	                    			 #min_dFF = min(dFF,na.rm=TRUE),
																	                    			 peak_key = paste(chemical_condition,trackROI_key,bin,sep="-"))


															    tmp_time_of_Glu_max<-unique(selected_Glu_data$time_of_max)
															    tmp_Glu_max_dFF<- unique(selected_Glu_data$max_dFF)
															    tmp_peak_key<-unique(selected_Glu_data$peak_key)

															    selected_Ca_data <- Ca_trace %>% ungroup() %>%
																							  dplyr::filter(chemical_condition == current_chem_condition, absoluteTime >= bin_start & absoluteTime <= bin_end) %>%
																							  group_by(chemical_condition,trackROI_key) %>%
																	                    	  mutate(window_start = first(absoluteTime),
																	                    	         window_end = last(absoluteTime),
																	                    			 bin = bin_number,
																	                    			 #max_dFF_JF = max(dFF,na.rm=TRUE),
																									 #peak_color = ifelse(max_dFF_JF >= hline_val, "blue","white"),
																	                    			 time_of_max = tmp_time_of_Glu_max,
																	                    			 max_dFF = tmp_Glu_max_dFF,
																	                    			 #min_dFF = min(dFF,na.rm=TRUE),
																	                    			 peak_key = tmp_peak_key,
																	                    			 )


		 														
														        selected_data <- bind_rows(selected_Glu_data,selected_Ca_data)


		 														selected_intervals[[j]] <- selected_data

			                         }

			                         
			                         trace_binned<- bind_rows(selected_intervals)

			                         trace_binned <- trace_binned %>% group_by(sensor,chemical_condition, trackROI_key, bin,peak_key) %>% mutate(normTime = absoluteTime - time_of_max)
			                         df_traces_annotated[[i]]<-trace_binned

			                      
			                     }
df_traces_annotated<- bind_rows(df_traces_annotated)
df_traces_annotated
}



#Glu_trace_binned <- bind_rows(selected_intervals) %>% group_by(chemical_condition,trackROI_key)
			                         #get_bins<- Glu_trace_binned %>% ungroup() %>% group_by(chemical_condition,trackROI_key) %>%  select(absoluteTime,bin, window_start,window_end,time_of_max,max_dFF,peak_key)
			                         
			                         #Ca_trace_binned<- suppressMessages(left_join(Ca_trace,get_bins,by=c("chemical_condition","trackROI_key","absoluteTime") )) #%>% dplyr::filter(!is.na(peak_key)) )

			                         

