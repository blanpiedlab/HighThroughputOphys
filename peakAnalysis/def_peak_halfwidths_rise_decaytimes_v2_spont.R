#def_peak_halfwidths_rise_decaytimes_v1.R
library(tictoc)


# in the summary dataframe, decay_hw - rise_hw will establish $halfwidth


get_peak_halfwidths_rise_decaytimes_v1<- function(df,groupers,remove_cols=drops) {

		
						### verifying that we only look at peaks that can be fit. 

						print("Generating the list of ROI_keys with taus.")

						lighter_df <- df[,!(names(df) %in% remove_cols)]
                        taus<-suppressMessages(lighter_df %>% group_by_at(groupers) %>%
                                     dplyr::filter(peakID != "NotPeak") %>%
                                     unnest(tidied) %>%
                                     spread(term,estimate) %>%
                                     ungroup() %>%
                                     group_by_at(groupers) %>%
                                     summarise( tau_decay = unique( tau_decay[!is.na(tau_decay)] ),                                                                              #tau_decays that are not NA, one for each peak. 
                                                A = unique( A[!is.na(A)] ) )
                                     )

                        ROIs_with_taus <- unique(taus$ROI_key)
                        #print(ROIs_with_taus)



                        print("Generating the parameters to calculate times from.")
						#precleaning
						clean_df <- lighter_df %>% 	group_by_at(groupers) %>%
											dplyr::filter(ROI_key %in% ROIs_with_taus) %>%
											#dplyr::filter(absoluteTime >= firstStim & absoluteTime < (firstStim + 0.5)) %>%
											dplyr::filter(peakID != "NotPeak") %>%
											mutate(special_normTime = absoluteTime - min(absoluteTime),
													max_dFF = max(dFF,na.rm=TRUE),
													max_time = ifelse(!is.na(max_dFF), absoluteTime[which(dFF == max_dFF)], NA),
													min_time = ifelse(!is.na(max_time), max_time - 0.025, absoluteTime[which(special_normTime==0)]) ) %>%
											dplyr::filter(absoluteTime >= unique(min_time), absoluteTime <= max_time) %>%
											summarise(max_dFF = unique(max_dFF),
													max_time = unique(max_time),
													zero_dFF = dFF[which.min(dFF)], 
													zero_time = absoluteTime[which.min(dFF)],
													pct10_dFF = ifelse(!is.na(max_dFF), 0.1*(max_dFF-zero_dFF), NA),
													pct90_dFF = ifelse(!is.na(max_dFF), 0.9*(max_dFF-zero_dFF), NA),
													half_dFF = ifelse(!is.na(max_dFF), 0.5*(max_dFF-zero_dFF), NA))

													#zero_dFF = ifelse(!is.na(rise_edge_window), dFF[which.min(dFF[which(absoluteTime <= max_time & absoluteTime >= rise_edge_window)])], NA), 
													#zero_time = ifelse(!is.na(zero_dFF),absoluteTime[which.min(abs(dFF - zero_dFF))], NA),
													

						print("Generating rise time data.")
						riseTimes<- clean_df %>% group_by_at(groupers) %>%
											#dplyr::filter(absoluteTime >= firstStim & absoluteTime <= max_time ) %>%
											dplyr::filter(!is.na(half_dFF)) %>%
											summarise(rise_slope = (unique(max_dFF) - unique(zero_dFF))/(unique(max_time)-unique(zero_time) ),
														rise90_time = unique(unique(pct90_dFF)/rise_slope)+zero_time, 
														rise10_time = unique(unique(pct10_dFF)/rise_slope)+zero_time,
														rise_time = rise90_time - rise10_time,
														half_rise_time = (unique(half_dFF)/rise_slope)+zero_time) 

						print("Generating decay time data.")

						#decay_half_amplitude = 0.5*(A),
					    #                     decay_half_time = maxTime+(-(log( decay_half_amplitude/A )*tau_decay) )

						decayTimes<- left_join(clean_df,taus) %>% group_by_at(groupers) %>%
											dplyr::filter(!is.na(half_dFF), !is.na(max_dFF)) %>%
											#dplyr::filter(special_normTime >= special_normTime[which.min(abs(dFF - max_dFF))]) %>%
											#mutate(pctile = percent_rank(dFF)) %>%
											#dplyr::filter(dFF <= unique(pct90_dFF) & dFF >= unique(pct10_dFF) ) %>%
											summarise(#tau_decay = unique(tau_decay),
														#decay10_time = absoluteTime[which.min(abs(pct10_dFF - dFF) )],
														#decay90_time = absoluteTime[which.min(abs(pct90_dFF - dFF) )],
														decay10_time = max_time+(-(log( pct10_dFF/A )*tau_decay) ),
														decay90_time = max_time+(-(log( pct90_dFF/A )*tau_decay) ),
														 
														decay_time = decay10_time - decay90_time,

														half_decay_time = max_time+(-(log( half_dFF/A )*tau_decay) ),
														
														#half_decay_time = absoluteTime[which.min(abs(half_dFF-dFF))],
														
														pct90_decay_dFF = pct90_dFF,#dFF[which(absoluteTime == decay90_time)],
														pct10_decay_dFF = pct10_dFF,#dFF[which(absoluteTime == decay10_time)], 
														half_decay_dFF = half_dFF)
														#dFF[which(absoluteTime == half_decay_time)])

						gc()
						join_df_a <- suppressMessages(left_join(clean_df, riseTimes))
						gc()
						join_df_all <- suppressMessages(left_join(join_df_a, decayTimes))

						rm(join_df_a, decayTimes, riseTimes)
						gc()
						print("outputting data frame summary.")

						output_df <- join_df_all %>% group_by_at(groupers) %>%
											 summarise(rise10_time = unique(rise10_time),
											 			rise90_time = unique(rise90_time),
											 			half_rise_time = unique(half_rise_time),
											 			rise_time = unique(rise_time),

											 			rise90_dFF = unique(pct90_dFF),
											 			rise10_dFF = unique(pct10_dFF),
											 			half_rise_dFF = unique(half_dFF),

											 			decay10_time = unique(decay10_time),
											 			decay90_time = unique(decay90_time),
											 			half_decay_time = unique(half_decay_time),
											 			decay_time = unique(decay_time),

											 			#tau_decay = unique(tau_decay),

											 			decay90_dFF = unique(pct90_decay_dFF),
											 			decay10_dFF = unique(pct10_decay_dFF),
											 			half_decay_dFF = unique(half_decay_dFF),

											 			halfwidth = half_decay_time-half_rise_time)
						rm(join_df_all)
						gc()

						output_df


}




#next steps in the pipeline
#which is better? half widths based on A or halfwidths based on true max amplitude? 
