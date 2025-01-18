

#subDir <- "peakStats_activity_vs_marker_v1"   


source(paste0(path,"peakAnalysis/def_interSpike_v1.R") )

source(paste0(path,"peakAnalysis/def_peak_halfwidths_rise_decaytimes_v1.R") )




#,plotBy,levels, color_override =NULL

peakStats<- function(df,groupers,subgroup,list_stimuli=stimList,keys, remove_cols=drops){
                        
                        

       # print("WARNING: INTERSPIKES ARE NOT BEING CALCULATED IN THIS ITERATION")
                        #### generate interSpikes
                        tic()
                        print("Generating interSpikes - Delta t.")
                        interSpike_df<- def_interSpike(df=df, groupers=groupers,keys=keys, list_stims = list_stimuli)
                        print("Finished generating interSpikes.")
                        toc()

                        ### generate peak amplitudes
                        tic()
                        print("Generating amplitudes.")

                        onlyCleanPeaks<- suppressMessages(interSpike_df %>% group_by_at(groupers) %>% dplyr::filter(peakID != "NotPeak") %>%
                                                            summarise(amplitude = max(dFF),
                                                                        interSpike = unique(interSpike)
                                                                        )
                                                        )

                        print("Finished generating amplitudes.")
                        toc()
                        

                        # #### generate taus
                        tic()
                        print("Generating tau_decays.")
                        lighter_df <- df[,!(names(df) %in% remove_cols)]
                        taus<-suppressMessages(lighter_df %>% group_by_at(groupers) %>%
                                     dplyr::filter(peakID != "NotPeak") %>%
                                     unnest(tidied) %>%
                                     spread(term,estimate) %>%
                                     ungroup() %>%
                                     group_by_at(groupers) %>%
                                     summarise( tau_decay = unique( tau_decay[!is.na(tau_decay)] ),                                                                              #tau_decays that are not NA, one for each peak. 
                                                A = unique( A[!is.na(A)] )                                                                              #tau_decays that are not NA, one for each peak. 
                                                )
                                     )

                        print("Finished generating tau_decays.")
                        toc()


                        print("Generating risetimes, decaytimes, and halfwidths.")


                        ROIs_with_taus <- unique(taus$ROI_key)




                        print("Generating the parameters to calculate times from.")
                        #precleaning
                        clean_df <- suppressMessages(lighter_df %>%  group_by_at(groupers) %>%
                                            dplyr::filter(ROI_key %in% ROIs_with_taus) %>%
                                            dplyr::filter(absoluteTime >= firstStim & absoluteTime < (firstStim + 0.5)) %>%
                                            dplyr::filter(peakID != "NotPeak") %>%
                                            mutate(special_normTime = absoluteTime - min(absoluteTime),
                                                    max_dFF = max(dFF,na.rm=TRUE),
                                                    max_time = ifelse(!is.na(max_dFF), absoluteTime[which(dFF == max_dFF)], NA),
                                                    min_time = ifelse(!is.na(max_time), absoluteTime[which(special_normTime==0)],NA) ) %>% ##max_time - 0.025 
                                             dplyr::filter(absoluteTime >= unique(min_time), absoluteTime <= max_time) %>%
                                            summarise(max_dFF = unique(max_dFF),
                                                    max_time = unique(max_time),
                                                    zero_dFF = dFF[which.min(dFF)], 
                                                    zero_time = absoluteTime[which.min(dFF)],
                                                    pct10_dFF = ifelse(!is.na(max_dFF), 0.1*(max_dFF), NA),
                                                    pct90_dFF = ifelse(!is.na(max_dFF), 0.9*(max_dFF), NA),
                                                    half_dFF = ifelse(!is.na(max_dFF), 0.5*(max_dFF), NA))
                                            )

                        print("Generating rise time data.")
                        riseTimes<- suppressMessages(clean_df %>% group_by_at(groupers) %>%
                                            #dplyr::filter(absoluteTime >= firstStim & absoluteTime <= max_time ) %>%
                                            dplyr::filter(!is.na(half_dFF)) %>%
                                             summarise(rise_slope = (unique(max_dFF) - unique(zero_dFF))/(unique(max_time)-unique(zero_time) ),
                                                        rise90_time = unique(unique(pct90_dFF)/rise_slope)+zero_time, 
                                                        rise10_time = unique(unique(pct10_dFF)/rise_slope)+zero_time,
                                                        rise_time = rise90_time - rise10_time,
                                                        half_rise_time = (unique(half_dFF)/rise_slope)+zero_time)
                                            ) 



                        print("Generating decay time data.")


                         print("Generating decay time data.")
                        decayTimes<- suppressMessages(left_join(clean_df,taus) %>% group_by_at(groupers) %>%
                                            #dplyr::filter(special_normTime >= special_normTime[which(dFF == max_dFF)]) %>%
                                            dplyr::filter(!is.na(half_dFF), !is.na(max_dFF)) %>%
                                            summarise(  decay10_time = max_time+(-(log( pct10_dFF/A )*tau_decay) ),
                                                        decay90_time = max_time+(-(log( pct90_dFF/A )*tau_decay) ),   
                                                        decay_time = decay10_time - decay90_time,

                                                        half_decay_time = max_time+(-(log( half_dFF/A )*tau_decay) ),
                                                        
                                                        pct90_decay_dFF = pct90_dFF,
                                                        pct10_decay_dFF = pct10_dFF, 
                                                        half_decay_dFF = half_dFF)
                                            )
                                                        

                        print("outputting data frame summary.")

                          gc()
                        join_df_a <- suppressMessages(left_join(onlyCleanPeaks, taus))
                        join_df_b <- suppressMessages(left_join(riseTimes, decayTimes))
                        
                        
                        merge_data_all <- suppressMessages(left_join(join_df_a,join_df_b) %>% group_by_at(groupers) %>%
                                             summarise(amplitude = unique(amplitude),
                                                        tau_decay = unique(tau_decay),
                                                        interSpike = unique(interSpike),
                                                        
                                                        rise10_time = unique(rise10_time),
                                                        rise90_time = unique(rise90_time),,
                                                        half_rise_time = unique(half_rise_time),

                                                        rise10_dFF = unique(pct10_decay_dFF),
                                                        rise90_dFF = unique(pct90_decay_dFF),
                                                        half_rise_dFF = unique(half_decay_dFF),

                                                        rise_time = unique(rise_time),

                                                        decay90_time = unique(decay90_time),
                                                        decay10_time = unique(decay10_time),
                                                        half_decay_time = unique(half_decay_time),

                                                        decay90_dFF = unique(pct90_decay_dFF),
                                                        decay10_dFF = unique(pct10_decay_dFF),
                                                        half_decay_dFF = unique(half_decay_dFF),

                                                        decay_time = unique(decay_time),

                                                        
                                                        
                                                        halfwidth = half_decay_time-half_rise_time) )
                                                
                        

                        rm(join_df_a,join_df_b,riseTimes,decayTimes,clean_df,onlyCleanPeaks,taus,lighter_df,ROIs_with_taus)

                        

                        tic()
                        print("Cleaning up peakStats for export.")
                        ### clean up peakStats
                        peakStats<- suppressMessages(merge_data_all %>% group_by_at(groupers) %>% 
                                                                summarise(
                                                                            tau_decay_ms = tau_decay*1000,
                                                                            interSpike_ms = interSpike*1000,
                                                                            amplitude = amplitude,
                                                                            
                                                                            rise10_time = rise10_time,
                                                                            rise10_dFF = rise10_dFF,
                                                                            
                                                                            rise90_time = rise90_time,
                                                                            rise90_dFF = rise90_dFF,
                                                                            
                                                                            half_rise_time = half_rise_time,
                                                                            half_rise_dFF = half_rise_dFF,

                                                                            t_rise = rise_time*1000,
                                                                            
                                                                            decay10_time = decay10_time,
                                                                            decay10_dFF = decay10_dFF,

                                                                            decay90_time = decay90_time,
                                                                            decay90_dFF = decay90_dFF,

                                                                            half_decay_time = half_decay_time,
                                                                            half_decay_dFF = half_decay_dFF,

                                                                            t_decay = decay_time*1000,
                                                                            t_half = halfwidth*1000) %>%
                                                                ungroup() %>%
                                                                group_by(Ca,segmentation) %>% 
                                                                 dplyr::filter( tau_decay_ms <= 200,
                                                                                tau_decay_ms >= 10, 
                                                                                ### DEPRECATED
                                                                                tau_decay_ms <= quantile(tau_decay_ms, 0.985, na.rm=TRUE), 
                                                                                tau_decay_ms >= quantile(tau_decay_ms, 0.015, na.rm=TRUE),
                                                                                
                                                                                interSpike_ms < 150,
                                                                                

                                                                                                        ### DEPRECATED
                                                                                                        #tau_decay_ms < quantile(tau_decay_ms, 0.985, na.rm=TRUE), 
                                                                                                        #tau_decay_ms > quantile(tau_decay_ms, 0.015, na.rm=TRUE),
                                                                                                        #amplitude < quantile(amplitude, 0.985, na.rm=TRUE), 
                                                                                                        #amplitude > quantile(amplitude, 0.015, na.rm=TRUE),
                                                                                                        #interSpike_ms < 210,

                                                                                amplitude >= 0.05,
                                                                                
                                                                                t_half > 0,
                                                                                t_rise > 0,
                                                                                t_decay > 0) 
                                                        )
                        print("Finished cleaning up peakStats for export.")
                        toc()


                        rm(merge_data_all)


                        
                   

peakStats
}
