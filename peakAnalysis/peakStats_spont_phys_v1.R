

#subDir <- "peakStats_activity_vs_marker_v1"   


#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/def_stim_vlines.R")
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/def_interSpike_v1.R")

#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/def_peak_halfwidths_rise_decaytimes_v1.R") 




#,plotBy,levels, color_override =NULL

peakStats<- function(df,groupers,remove_cols=drops){
                        
                            
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
                     


                        print("Generating the parameters to calculate times from.")
                        #precleaning
                        clean_df <- suppressMessages(lighter_df %>%  group_by_at(groupers) %>%
                                            dplyr::filter(ROI_key %in% ROIs_with_taus) %>%
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
                                            )
                        print("Extracting maxima.")
                        onlyMax <- suppressMessages(clean_df %>% group_by_at(groupers) %>% summarise(amplitude = unique(max_dFF)) )
                                                    

                        print("Generating rise time data.")
                        riseTimes<- suppressMessages(clean_df %>% group_by_at(groupers) %>%
                                            dplyr::filter(!is.na(half_dFF)) %>%
                                            summarise(rise_slope = (unique(max_dFF) - unique(zero_dFF))/(unique(max_time)-unique(zero_time) ),
                                                        rise90_time = unique(unique(pct90_dFF)/rise_slope)+zero_time, 
                                                        rise10_time = unique(unique(pct10_dFF)/rise_slope)+zero_time,
                                                        rise_time = rise90_time - rise10_time,
                                                        half_rise_time = (unique(half_dFF)/rise_slope)+zero_time)
                                            ) 



                        print("Generating decay time data.")
                        decayTimes<- suppressMessages(left_join(clean_df,taus) %>% group_by_at(groupers) %>%
                                            dplyr::filter(!is.na(half_dFF), !is.na(max_dFF)) %>%
                                            summarise(  decay10_time = max_time+(-(log( pct10_dFF/A )*tau_decay) ),
                                                        decay90_time = max_time+(-(log( pct90_dFF/A )*tau_decay) ),   
                                                        decay_time = decay10_time - decay90_time,

                                                        half_decay_time = max_time+(-(log( half_dFF/A )*tau_decay) ),
                                                        
                                                        pct90_decay_dFF = pct90_dFF,
                                                        pct10_decay_dFF = pct10_dFF, 
                                                        half_decay_dFF = half_dFF)
                                            )
                                                        
                        gc()
                        join_df_a <- suppressMessages(left_join(onlyMax, taus))
                        join_df_b <- suppressMessages(left_join(riseTimes, decayTimes))
                        
                        
                        merge_data_all <- suppressMessages(left_join(join_df_a,join_df_b) %>% group_by_at(groupers) %>%
                                             summarise(amplitude = unique(amplitude),
                                                        tau_decay = unique(tau_decay),
                                                        #rise10_time = unique(rise10_time),
                                                        #rise90_time = unique(rise90_time),
                                                        half_rise_time = unique(half_rise_time),
                                                        rise_time = unique(rise_time),

                                                        #rise90_dFF = unique(pct90_dFF),
                                                        #rise10_dFF = unique(pct10_dFF),
                                                        #half_rise_dFF = unique(half_dFF),

                                                        #decay10_time = unique(decay10_time),
                                                        #decay90_time = unique(decay90_time),
                                                        half_decay_time = unique(half_decay_time),
                                                        decay_time = unique(decay_time),

                                                        #tau_decay = unique(tau_decay),

                                                        #decay90_dFF = unique(pct90_decay_dFF),
                                                        #decay10_dFF = unique(pct10_decay_dFF),
                                                        #half_decay_dFF = unique(half_decay_dFF),

                                                        halfwidth = half_decay_time-half_rise_time) )
                                                
                        rm(join_df_a,join_df_b,riseTimes,decayTimes,clean_df,onlyMax,taus,lighter_df,ROIs_with_taus)

                        gc()
   
                        
                        tic()
                        print("Cleaning up peakStats for export.")
                        ### clean up peakStats
                        peakStats<- suppressMessages(merge_data_all %>% group_by_at(groupers) %>% 
                                                                summarise(tau_decay_ms = tau_decay*1000,
                                                                            amplitude = amplitude,
                                                                            t_rise = rise_time*1000,
                                                                            t_decay = decay_time*1000,
                                                                            t_half = halfwidth*1000) %>%
                                                                ungroup() %>%
                                                                group_by(Ca,segmentation) %>%
                                                                #group_by_at(subgroup) %>% 
                                                                 dplyr::filter( tau_decay_ms <= 150,
                                                                                tau_decay_ms >= 10, 
                                                                                
                                                                                #tau_decay_ms < quantile(tau_decay_ms, 0.985, na.rm=TRUE), 
                                                                                #tau_decay_ms > quantile(tau_decay_ms, 0.015, na.rm=TRUE),
                                                                                
                                                                                #amplitude < quantile(amplitude, 0.985, na.rm=TRUE), 
                                                                                #amplitude > quantile(amplitude, 0.015, na.rm=TRUE),
                                                                                amplitude >= 0.1,
                                                                                t_half > 0) 
                                                        )
                        print("Finished cleaning up peakStats for export.")
                        toc()


                        rm(merge_data_all)


                        
                   

peakStats
}









