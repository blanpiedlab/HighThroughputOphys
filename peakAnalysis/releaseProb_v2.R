##BasalRP_stats_v2.R
## Written by Samuel T. Barlow
## 10.19.22

# BasalRP_stats_v2.R is a module designed to capitalize on the richness of data that can be extracted from marker-based segmentation experiments. 
# For SfN 2022, my goal is to use marker-based segmentation to extract the following properties from the data-stream: 
# 1. Calculate pseudo-release probability
# 2. Calculate quantal size (will require spont data)
# 3. Calculate coefficient of variation per ROI
# 4. Averaged traces 
# 5. generate a per-ROI release prob, quantal size map visualization a la PPR_plotter.R  






#library(M3C)



subDir <- "track_singleAP_phys_v10_reviseUMAPs_again"   


source(paste0(path,"sub-functions/setFigureDirectory.R") ) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"sub-functions/def_stim_vlines.R"))
source(paste0(path,"peakAnalysis/png_plotter_v2_basalRP.R"))
#source(paste0(path,"peakAnalysis/png_plotter_v3_ROIoverlay.R"))
#source(paste0(path,"peakAnalysis/png_plotter_v4.R"))
source(paste0(path,"peakAnalysis/find_replicates.R"))
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))




library(ggforce)
library(scales)
library(nls.multstart)

library(umap)
library(dbscan)
library(RColorBrewer)


releaseProb_stats<- function(peak_stat_df = peak_stats, df = df_traces, spont_avgs = spont_avg_df, MV_analysis = MV_output, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops,ROIs_to_save = ROIs_to_save){




                        #### establish variables that will be used for generating datasets


                        vid_keys_to_save = unique(str_extract(ROIs_to_save, "dish\\d\\d-plate\\d\\d-region\\d\\d") ) 
                        color_switch = !is.null(color_override)
                        plot_height = 8
                        interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE) 
                        


                        remove_neurons = c("dish05-plate03-region02", "dish05-plate06-region01", "dish05-plate06-region02")#,"dish05-plate01-region01")
                        #### clean up dataframe by merging with peakStats

                        peak_stat_df <- peak_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID), 
                                                                Ca_mM = case_when(Ca == "0pt5Ca" ~ 0.5,
                                                                                  Ca == "1Ca" ~ 1,
                                                                                  Ca == "2Ca" ~ 2,
                                                                                  Ca == "4Ca" ~ 4),
                                                                imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>%
                                                  dplyr::filter(!imaging_region %in% remove_neurons)
                        
                        light_df <- df[,!(names(df) %in% remove_col)]

                       
                        ROIs_to_remove <- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are finite
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%
                                                 dplyr::filter( is.na(dFF))
                        ROIs_to_remove_timeskip<- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are do not have a frame skip in time of interest
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%            
                                                 summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
                                                 dplyr::filter(maxInterframe > 0.18) %>%
                                                 select(ROI_key,maxInterframe) 

#print(paste0("These are the bad recordings: ", unique(ROIs_to_remove$ROI_key)))
#print(paste0("These are the bad recordings with a frame jump during post-stimulus: ", unique(ROIs_to_remove_timeskip$ROI_key)))

                       
                        # print("These should be the borked ROIs we found (drift correction failures?)")
                        # print(unique(ROIs_to_remove$ROI_key))
                        # check<- ROIs_to_remove
                        ROIs_to_remove <- unique(ROIs_to_remove$ROI_key)     
                        ROIs_to_remove_timeskip<-unique(ROIs_to_remove_timeskip$ROI_key)                                             

                        clean_df <- light_df %>%  dplyr::filter(!ROI_key %in% ROIs_to_remove, !ROI_key %in% ROIs_to_remove_timeskip) %>% 
                                                  mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>%
                                                  dplyr::filter(!imaging_region %in% remove_neurons)
                        
                        replicate_groupers= c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","vid_key","trackROI_key")
                        colname_toFind = "exposeNum"
                        df_mutated<- find_replicates(df = clean_df, groupers=replicate_groupers, colname_toFind=colname_toFind)
                        clean_df<-df_mutated
                        print("These were the replicates we found after cleaning")
                        print(unique(clean_df$replicates))
                        rm(df_mutated)


                        clean_peaks<- suppressMessages(left_join(clean_df,peak_stat_df) %>% dplyr::filter(windowedPeakID != "NotPeak") ) 
                        timezero_df<- suppressMessages(clean_peaks %>%  group_by_at(groupers) %>% summarise(time_zero = unique(firstStim[!is.na(firstStim)] ) ) )
                        clean_peaks_normTime<- suppressMessages(left_join(clean_peaks,timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))
                        

                        
                        new_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
                        
                        clean_peaks_normTime$Ca<- as.factor(as.character( (clean_peaks_normTime$Ca) ) )
                        peak_stat_df$Ca<- as.factor(as.character( (peak_stat_df$Ca) ) )
                        

                        levels(peak_stat_df$Ca) <- new_levels
                        levels(clean_peaks_normTime$Ca) <- new_levels

                        
                      

                        #### generate ROI count from raw data

                        getROIs<- clean_df %>% group_by(vid_key,Ca, Ca_mM) %>% summarise(ROI_count = length(unique(ROI_key)))
                        get_trackedROIs <- clean_df %>% group_by(Ca) %>% mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key) ) %>% 
                                                                        summarise(imaging_region_count = length(unique(imaging_region)),
                                                                                    trackedROI_count  = length(unique(trackROI_key)))                        

                        area_imaging_region_um = 25.6 * 25.6
                        check_ROIs_per_vid <- suppressMessages(getROIs %>% group_by(Ca) %>% mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>% 
                                                                                                          summarise(total_imaging_regions = length(unique(imaging_region)),
                                                                                                                        total_vids = length(unique(vid_key)),
                                                                                                                        sum_ROIs = sum(ROI_count),
                                                                                                                        ROIs_per_vid = sum_ROIs/total_vids,
                                                                                                                        mean_ROI_count_per_vid = mean(ROI_count),
                                                                                                                        sd_ROI_count = sd(ROI_count),
                                                                                                                        se_ROI_count = sd_ROI_count/sqrt(sum_ROIs),
                                                                                                                        as_area_mean_ROI = mean_ROI_count_per_vid/area_imaging_region_um,
                                                                                                                        as_area_sd_ROI = sd_ROI_count/area_imaging_region_um ) )

                        #print(getROIs)






    ROIs_with_peaks<- suppressMessages(clean_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5), !ROI_key %in% ROIs_to_remove) %>%  
                                                                group_by(Ca, Ca_mM, vid_key,trackROI_key,ROI_key,ROINumber,exposeNum,replicates) %>%  
                                                                summarise(peak_detected = length(unique(peakID))-1,
                                                                                peak_positive = case_when(peak_detected > 0 ~ 1,
                                                                                                          peak_detected == 0 ~ 0) ) #%>%
                                                                #dplyr::filter(peak_detected <= 2) # get rid of those with multiple peaks
                                                                )

                    droppers = c("sensor","marker","segmentation","TTL_start","dish","plate","region","peakID","protocol",'exposeNum', "windowedPeakID","Ca")
                    cleanup_peak_stat_data_point <- peak_stat_df %>% dplyr::filter(ROI_key == "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005", interSpike_ms < 100)
                    cleanup_peak_stat_df <- peak_stat_df %>% dplyr::filter(ROI_key != "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005" )
                    cleanup_merge <- bind_rows(cleanup_peak_stat_df, cleanup_peak_stat_data_point)
                    light_peak_stat_df <- cleanup_merge[,!(names(cleanup_merge) %in% droppers)]

                    #spont_avg_data for quantal size calc
                    light_spont_avg <- spont_avgs %>% mutate(Ca_mM = case_when(Ca == "\"0.5 mM \" * Ca^\"2+\"" ~ 0.5,
                                                                                    Ca == "\"1 mM \" * Ca^\"2+\"" ~  1,
                                                                                    Ca == "\"2 mM \" * Ca^\"2+\"" ~ 2,
                                                                                    Ca == "\"4 mM \" * Ca^\"2+\"" ~ 4),
                                                            spont_amplitude = mean_amplitude ) %>%
                                                    ungroup() %>%
                                                    select(Ca_mM,spont_amplitude)
                    #print("Here's a glimpse of the spont_avg_df")
                    #print(head(light_spont_avg))

                    light_peak_stat_df <- suppressMessages(left_join(light_peak_stat_df,light_spont_avg) %>% 
                                                                            mutate(quanta_calc = amplitude/spont_amplitude, 
                                                                                    Ca_expr = case_when(Ca_mM == 0.5 ~ "0pt5Ca",
                                                                                                    Ca_mM == 1 ~ "1Ca",
                                                                                                    Ca_mM ==  2 ~ "2Ca",
                                                                                                    Ca_mM == 4 ~ "4Ca")) )
                    print("Here's a glimpse of the quanta_calc")
                    print(unique(light_peak_stat_df$quanta_calc))
                    light_peak_stat_df$Ca_expr <- factor(light_peak_stat_df$Ca_expr, levels = c("0pt5Ca","1Ca","2Ca","4Ca"))
                    levels(light_peak_stat_df$Ca_expr) <- new_levels

                                                            
                    

                    print("These are the columns that exist in peak_stat_df")
                    #print(names(peak_stat_df))
                    print(names(light_peak_stat_df))

                    print(head(light_peak_stat_df))



####### figure out new mean_variance analysis


                    # estimate_max_signal_df<- clean_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.15)) %>%
                    #                                      group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>% 
                    #                                      summarise(amplitude = max(dFF, na.rm=TRUE))

                    # mean_noise_df <- clean_df %>% dplyr::filter(peakID == "NotPeak") %>%
                    #                                      group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>% 
                    #                                      summarise(median_noise = median(dFF, na.rm=TRUE))

                    # estimates_per_trial_df <- left_join(estimate_max_signal_df,mean_noise_df) %>% mutate(normalize_amplitude = amplitude - median_noise)
                    # estimates_per_trial_df <- left_join(estimates_per_trial_df, light_spont_avg) %>% mutate(corrected_amplitude = case_when(normalize_amplitude < spont_amplitude ~ 0,
                    #                                                                                                                          normalize_amplitude >= spont_amplitude ~ normalize_amplitude),
                    #                                                                                         quanta_calc = corrected_amplitude/spont_amplitude)

                    # MV_per_synapse_df <- estimates_per_trial_df %>% group_by(Ca, Ca_mM, trackROI_key) %>% 
                    #                                                 summarise(n_obs = n(),
                    #                                                           mean_amplitude = mean(corrected_amplitude,na.rm=TRUE),
                    #                                                           var_amplitude = var(corrected_amplitude,na.rm=TRUE),
                    #                                                           mean_quanta = mean(quanta_calc, na.rm=TRUE), 
                    #                                                           var_quanta = var(quanta_calc,na.rm=TRUE),
                    #                                                           sd_amplitude = sd(corrected_amplitude,na.rm=TRUE),
                    #                                                           se_amplitude = sd_amplitude/sqrt(n_obs),
                    #                                                           CV_amplitude = sd_amplitude/mean_amplitude

                    #                                                           )
                                                                     

                    # keep_ROIs_MV <- MV_per_synapse_df %>% dplyr::filter( Ca == "4Ca")
                    # keep_ROIs_MV_list<- unique(keep_ROIs_MV$trackROI_key)
                    # MV_per_synapse_df <- MV_per_synapse_df %>% dplyr::filter(trackROI_key %in% keep_ROIs_MV_list)




















###### CALCULATE RELEASE PROBABILITY AND OTHER PARAMS PER ROI #######
                    

                    

                
                         
                    join_df <- left_join(ROIs_with_peaks, light_peak_stat_df)

                    print(head(join_df))
                    releaseProb_calc <- suppressMessages( join_df %>% 
                                                         ungroup() %>%
                                                         group_by(Ca,Ca_mM,vid_key,trackROI_key,ROINumber) %>%
                                                         summarise(n_peaks = sum(peak_detected),
                                                                    releaseProbability_perROI = ifelse(n_peaks/first(unique(replicates))<=1, n_peaks/first(unique(replicates)), 1), #,#n_peaks/first(unique(replicates)),#
                                                                    activity = case_when(releaseProbability_perROI ==  0 ~ "inactive",
                                                                                         releaseProbability_perROI >= 0.75  ~ "very active",
                                                                                         releaseProbability_perROI >= 0.50 ~ "active",
                                                                                         releaseProbability_perROI >= 0.25 ~ "somewhat active",
                                                                                         releaseProbability_perROI > 0 ~ "rarely active"#,
                                                                                         #releaseProbability_perROI > 1.0 ~ "multivesicular"
                                                                                         ),
                                                                    mean_quanta = mean(quanta_calc, na.rm=TRUE), 
                                                                    mean_amplitude = mean(amplitude,na.rm=TRUE),
                                                                    sd_amplitude = sd(amplitude,na.rm=TRUE),
                                                                    se_amplitude = sd_amplitude/sqrt(n_peaks),
                                                                    var_amplitude = var(amplitude,na.rm=TRUE),
                                                                    var_quanta = var(quanta_calc,na.rm=TRUE),
                                                                    CV_amplitude = sd_amplitude/mean_amplitude,

                                                                    
                                                                   
                                                                    mean_dt = mean(interSpike_ms,na.rm=TRUE),
                                                                    sd_dt = sd(interSpike_ms,na.rm=TRUE),
                                                                    se_dt = sd_dt/sqrt(n_peaks),
                                                                    var_dt = var(interSpike_ms, na.rm=TRUE),

                                                                    CV_dt = sd_dt/mean_dt,

                                                                    mean_tau_decay = mean(tau_decay_ms,na.rm=TRUE),
                                                                    sd_tau_decay = sd(tau_decay_ms,na.rm=TRUE),
                                                                    se_tau_decay = sd_tau_decay/sqrt(n_peaks),
                                                                    var_tau_decay = var(tau_decay_ms, na.rm=TRUE),

                                                                    CV_tau_decay = sd_tau_decay/mean_tau_decay,

                                                                    mean_t_rise = mean(t_rise,na.rm=TRUE),
                                                                    sd_t_rise = sd(t_rise,na.rm=TRUE),
                                                                    se_t_rise = sd_t_rise/sqrt(n_peaks),

                                                                    CV_t_rise = sd_t_rise/mean_t_rise,

                                                                    mean_t_decay = mean(t_decay,na.rm=TRUE),
                                                                    sd_t_decay = sd(t_decay,na.rm=TRUE),
                                                                    se_t_decay = sd_t_decay/sqrt(n_peaks),

                                                                    CV_t_decay = sd_t_decay/mean_t_decay,
                                                                    
                                                                    mean_t_half = mean(t_half,na.rm=TRUE),
                                                                    sd_t_half = sd(t_half,na.rm=TRUE),
                                                                    se_t_half = sd_t_half/sqrt(n_peaks),

                                                                    CV_t_half = sd_t_half/mean_t_half ) %>%
                                                         mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key))
                                                                    )



                    #### add in data from MV_per_synapse_df
#                        releaseProb_calc <- left_join(MV_per_synapse_df, releaseProb_calc)

                        
                        plot_avgs <- releaseProb_calc %>% group_by(Ca,Ca_mM) %>% 
                                                            summarise(n_obs = n(),
                                                                    mean_P_iGlu = mean(releaseProbability_perROI,na.rm=TRUE),
                                                                    sd_P_iGlu = sd(releaseProbability_perROI,na.rm=TRUE),
                                                                    se_P_iGlu = sd_P_iGlu/sqrt(n_obs),

                                                                    median_mean_quanta = median(mean_quanta, na.rm=TRUE),

                                                                    mean_mean_quanta = mean(mean_quanta,na.rm=TRUE),
                                                                    sd_mean_quanta = sd(mean_quanta,na.rm=TRUE),
                                                                    se_mean_quanta = sd_mean_quanta/sqrt(n_obs))



print("Finished calculating tracked physiology parameters.")                                                                            
print(head(releaseProb_calc))


                                                                            
   my_medians <- light_peak_stat_df %>%
                                        group_by(Ca_expr) %>%
                                      summarize(median_amplitude = median(amplitude,na.rm=TRUE),
                                                median_tau = median(tau_decay_ms,na.rm=TRUE),
                                                median_interSpike =median(interSpike_ms,na.rm=TRUE),
                                                median_t_rise = median(t_rise,na.rm=TRUE),
                                                        median_t_half = median(t_half,na.rm=TRUE),
                                                        median_t_decay = median(t_decay,na.rm=TRUE),
                                                        median_quanta = median(quanta_calc, na.rm=TRUE) ) #,

print("finished calculating median values.")


#### fix corr_stats to capture all stats



active_4Ca_list<- releaseProb_calc %>% dplyr::filter(Ca_mM == 4, releaseProbability_perROI > 0, !is.na(var_amplitude) )
active_4Ca_list<- unique(active_4Ca_list$trackROI_key)

releaseProb_calc <- releaseProb_calc %>% dplyr::filter( trackROI_key %in% active_4Ca_list)


#Generate RRP vs. PiGlu at 0.5 mM
#### mean_quanta
get_RRP_05Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_1Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_2Ca <- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_4Ca <- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)

names(get_RRP_05Ca)[names(get_RRP_05Ca) == 'mean_quanta'] <- 'mean_quanta_05Ca'
names(get_RRP_1Ca)[names(get_RRP_1Ca) == 'mean_quanta'] <- 'mean_quanta_1Ca'
names(get_RRP_2Ca)[names(get_RRP_2Ca) == 'mean_quanta'] <- 'mean_quanta_2Ca'
names(get_RRP_4Ca)[names(get_RRP_4Ca) == 'mean_quanta'] <- 'mean_quanta_4Ca'



#### var_quanta
get_var_quanta_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,var_quanta)
get_var_quanta_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,var_quanta)
get_var_quanta_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,var_quanta)
get_var_quanta_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,var_quanta)

names(get_var_quanta_05Ca)[names(get_var_quanta_05Ca) == 'var_quanta'] <- 'var_quanta_05Ca'
names(get_var_quanta_1Ca)[names(get_var_quanta_1Ca) == 'var_quanta'] <- 'var_quanta_1Ca'
names(get_var_quanta_2Ca)[names(get_var_quanta_2Ca) == 'var_quanta'] <- 'var_quanta_2Ca'
names(get_var_quanta_4Ca)[names(get_var_quanta_4Ca) == 'var_quanta'] <- 'var_quanta_4Ca'



###mean amplitude
get_dFF_05Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_1Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_2Ca <- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_4Ca <- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)

names(get_dFF_05Ca)[names(get_dFF_05Ca) == 'mean_amplitude'] <- 'mean_amplitude_05Ca'
names(get_dFF_1Ca)[names(get_dFF_1Ca) == 'mean_amplitude'] <- 'mean_amplitude_1Ca'
names(get_dFF_2Ca)[names(get_dFF_2Ca) == 'mean_amplitude'] <- 'mean_amplitude_2Ca'
names(get_dFF_4Ca)[names(get_dFF_4Ca) == 'mean_amplitude'] <- 'mean_amplitude_4Ca'


####CV_amplitude

get_CV_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)


names(get_CV_05Ca)[names(get_CV_05Ca) == 'CV_amplitude'] <- 'CV_amplitude_05Ca'
names(get_CV_1Ca)[names(get_CV_1Ca) == 'CV_amplitude'] <- 'CV_amplitude_1Ca'
names(get_CV_2Ca)[names(get_CV_2Ca) == 'CV_amplitude'] <- 'CV_amplitude_2Ca'
names(get_CV_4Ca)[names(get_CV_4Ca) == 'CV_amplitude'] <- 'CV_amplitude_4Ca'


#### var_amplitude
get_var_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)

names(get_var_05Ca)[names(get_var_05Ca) == 'var_amplitude'] <- 'var_amplitude_05Ca'
names(get_var_1Ca)[names(get_var_1Ca) == 'var_amplitude'] <- 'var_amplitude_1Ca'
names(get_var_2Ca)[names(get_var_2Ca) == 'var_amplitude'] <- 'var_amplitude_2Ca'
names(get_var_4Ca)[names(get_var_4Ca) == 'var_amplitude'] <- 'var_amplitude_4Ca'


#### PiGlu
get_PiGlu_05<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_1<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_2<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_4<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)

names(get_PiGlu_05)[names(get_PiGlu_05) == 'releaseProbability_perROI'] <- 'PiGlu_05Ca'
names(get_PiGlu_1)[names(get_PiGlu_1)   == 'releaseProbability_perROI'] <- 'PiGlu_1Ca'
names(get_PiGlu_2)[names(get_PiGlu_2)   == 'releaseProbability_perROI'] <- 'PiGlu_2Ca'
names(get_PiGlu_4)[names(get_PiGlu_4)   == 'releaseProbability_perROI'] <- 'PiGlu_4Ca'


#####mean_tau_decay
get_mean_tau_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)


names(get_mean_tau_05Ca)[names(get_mean_tau_05Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_05Ca'
names(get_mean_tau_1Ca)[names(get_mean_tau_1Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_1Ca'
names(get_mean_tau_2Ca)[names(get_mean_tau_2Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_2Ca'
names(get_mean_tau_4Ca)[names(get_mean_tau_4Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_4Ca'




#####CV_tau_decay
get_CV_tau_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)


names(get_CV_tau_05Ca)[names(get_CV_tau_05Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_05Ca'
names(get_CV_tau_1Ca)[names(get_CV_tau_1Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_1Ca'
names(get_CV_tau_2Ca)[names(get_CV_tau_2Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_2Ca'
names(get_CV_tau_4Ca)[names(get_CV_tau_4Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_4Ca'



####mean_dt

get_mean_dt_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)


names(get_mean_dt_05Ca)[names(get_mean_dt_05Ca) == 'mean_dt'] <- 'mean_dt_05Ca'
names(get_mean_dt_1Ca)[names(get_mean_dt_1Ca) == 'mean_dt'] <- 'mean_dt_1Ca'
names(get_mean_dt_2Ca)[names(get_mean_dt_2Ca) == 'mean_dt'] <- 'mean_dt_2Ca'
names(get_mean_dt_4Ca)[names(get_mean_dt_4Ca) == 'mean_dt'] <- 'mean_dt_4Ca'



#####CV_dt
get_CV_dt_05Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt5Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)


names(get_CV_dt_05Ca)[names(get_CV_dt_05Ca) == 'CV_dt'] <- 'CV_dt_05Ca'
names(get_CV_dt_1Ca)[names(get_CV_dt_1Ca) == 'CV_dt'] <- 'CV_dt_1Ca'
names(get_CV_dt_2Ca)[names(get_CV_dt_2Ca) == 'CV_dt'] <- 'CV_dt_2Ca'
names(get_CV_dt_4Ca)[names(get_CV_dt_4Ca) == 'CV_dt'] <- 'CV_dt_4Ca'





#### incorporate MV_params

MV_params<- MV_analysis$clean_N_sites

MV_synapse_list<- unique(MV_params$trackROI_key)

get_N_sites<- MV_params %>% dplyr::filter(Ca_mM == 4) %>% ungroup() %>% select(trackROI_key,N)

names(get_N_sites)[names(get_N_sites) == 'N'] <- 'N_sites'

get_Pr_05<- MV_params %>% dplyr::filter(Ca_mM == 0.5) %>% ungroup() %>% select(trackROI_key,P_rw)
get_Pr_1<- MV_params %>% dplyr::filter(Ca_mM == 1) %>% ungroup() %>% select(trackROI_key,P_rw)
get_Pr_2<- MV_params %>% dplyr::filter(Ca_mM == 2) %>% ungroup() %>% select(trackROI_key,P_rw)
get_Pr_4<- MV_params %>% dplyr::filter(Ca_mM == 4) %>% ungroup() %>% select(trackROI_key,P_rw)

names(get_Pr_05)[names(get_Pr_05) == 'P_rw'] <- 'Pr_05Ca'
names(get_Pr_1)[names(get_Pr_1) == 'P_rw'] <- 'Pr_1Ca'
names(get_Pr_2)[names(get_Pr_2) == 'P_rw'] <- 'Pr_2Ca'
names(get_Pr_4)[names(get_Pr_4) == 'P_rw'] <- 'Pr_4Ca'


get_RRP_ratio_05<- MV_params %>% dplyr::filter(Ca_mM == 0.5) %>% ungroup() %>% select(trackROI_key,RRP_ratio)
get_RRP_ratio_1<- MV_params %>% dplyr::filter(Ca_mM == 1) %>% ungroup() %>% select(trackROI_key,RRP_ratio)
get_RRP_ratio_2<- MV_params %>% dplyr::filter(Ca_mM == 2) %>% ungroup() %>% select(trackROI_key,RRP_ratio)
get_RRP_ratio_4<- MV_params %>% dplyr::filter(Ca_mM == 4) %>% ungroup() %>% select(trackROI_key,RRP_ratio)

names(get_RRP_ratio_05)[names(get_RRP_ratio_05) == 'RRP_ratio'] <- 'RRP_ratio_05Ca'
names(get_RRP_ratio_1)[names(get_RRP_ratio_1) == 'RRP_ratio'] <- 'RRP_ratio_1Ca'
names(get_RRP_ratio_2)[names(get_RRP_ratio_2) == 'RRP_ratio'] <- 'RRP_ratio_2Ca'
names(get_RRP_ratio_4)[names(get_RRP_ratio_4) == 'RRP_ratio'] <- 'RRP_ratio_4Ca'



get_Q<- MV_params %>% dplyr::filter(Ca_mM == 4) %>% ungroup() %>% select(trackROI_key,Q_w)

names(get_Q)[names(get_Q) == 'Q_w'] <- 'Q'







MV_list<- list(get_N_sites,
                get_Pr_05, get_Pr_1, get_Pr_2, get_Pr_4,
                get_RRP_ratio_05, get_RRP_ratio_1, get_RRP_ratio_2, get_RRP_ratio_4,
                get_Q
                )



#### 0.5 Ca often had zeros in the dataset not agnostic to signal definition. Omitted to avoid over-weighting insignificant properties
#get_dFF_05Ca, 
                 #get_dFF_1Ca, get_dFF_2Ca, get_dFF_4Ca, 
                 #get_var_05Ca, 
                 #get_var_1Ca,get_var_2Ca,get_var_4Ca,



stat_list<- list(get_RRP_05Ca, get_RRP_1Ca,get_RRP_2Ca,get_RRP_4Ca,
                 get_var_quanta_05Ca,
                 get_var_quanta_1Ca, get_var_quanta_2Ca, get_var_quanta_4Ca,
                 
                
                
                #get_mean_dt_05Ca, 
                #get_mean_dt_1Ca, get_mean_dt_2Ca, get_mean_dt_4Ca,
                
                #get_mean_tau_05Ca, 
                get_mean_tau_1Ca, get_mean_tau_2Ca, get_mean_tau_4Ca,
                
                get_PiGlu_05, get_PiGlu_1, get_PiGlu_2, get_PiGlu_4, 
                
                #get_CV_05Ca, 
                get_CV_1Ca, get_CV_2Ca, get_CV_4Ca, 
                 
                #get_CV_dt_05Ca, 
                get_CV_dt_1Ca, get_CV_dt_2Ca, get_CV_dt_4Ca,

                #get_CV_tau_05Ca,  
                get_CV_tau_1Ca, get_CV_tau_2Ca, get_CV_tau_4Ca
                )

corr_stats<- stat_list %>% reduce(full_join, by='trackROI_key') #%>% dplyr::filter(trackROI_key %in% active_4Ca_list)



#corr_stats<- corr_stats %>% #mutate(RRP_score = case_when( mean_quanta_4Ca >= quantile(mean_quanta_4Ca, 0.50, na.rm=TRUE) ~ "red", # ##mean_quanta_4Ca - mean_quanta_2Ca > 0
   #                                                                                mean_quanta_4Ca <= quantile(mean_quanta_4Ca, 0.50, na.rm=TRUE) ~ "blue"), #mean_quanta_4Ca/mean_quanta_2Ca ## #mean_quanta_4Ca - mean_quanta_2Ca <= 0
                                   #release_ratio = mean_quanta_4Ca / mean_quanta_2Ca
  #                                                                               ) %>%
 #                       dplyr::filter(!is.na(RRP_score))


MV_stats<- MV_list %>% reduce(full_join, by='trackROI_key')

corr_stats_w_MV<- left_join(MV_stats, corr_stats)







### UMAP code starting point


is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

standardize <- function(x){(x-min(x))/(max(x)-min(x))}


set.seed(333)

#get_UMAP is a function that accepts a list of stats to standardize and correlate for clustering in a UMAP


generate_UMAP<- function(stat_df,exclude_cols = exclude_cols){

            data_for_umap <- stat_df %>% select( -trackROI_key)

            data_for_umap[is.na(data_for_umap)] <- 0
            data_for_umap[is.nan(data_for_umap)] <- 0

            length_df_cols<- ncol(data_for_umap)
            
            standardized_umap<- data_for_umap %>% mutate_at(vars(-exclude_cols),  standardize) 



            synapse.umap <- umap(standardized_umap)
            df_xy<- as.data.frame(synapse.umap$layout)
            colnames(df_xy)[1]<- "x"
            colnames(df_xy)[2]<- "y"
            df_xy<- df_xy %>% mutate(index = row_number())
            df_data<- as.data.frame(synapse.umap$data)
            df_data<-df_data %>% mutate(index = row_number())
            df_annotated<-join(df_xy,df_data)


            df_annotated

}


cluster_UMAP<- function(UMAP_df, minPts = 10){


                    UMAP_xy_coords<- UMAP_df %>% select(x, y)
                    clustered_UMAP<- hdbscan(UMAP_xy_coords, minPts = minPts)
                    cluster_ID<- clustered_UMAP$cluster
                    UMAP_df_cluster_ID<- UMAP_df %>% bind_cols(cluster_ID = cluster_ID)
                    UMAP_df_cluster_ID$cluster_ID <- as.character(UMAP_df_cluster_ID$cluster_ID)
                    
                    UMAP_df_cluster_ID
 }


 
plot_UMAP<- function(UMAP_df, colour_var, plot_title, plot_filename, tag_var, global_output =TRUE, scale_bool = TRUE, colour_override =NULL){


                            x_min = min(UMAP_df$x)
                            x_max = max(UMAP_df$x)

                            y_min = min(UMAP_df$y)
                            y_max = max(UMAP_df$y)


                            if(scale_bool == TRUE) {
                                print("switching UMAP to a manual scale")

                            }
                            if(colour_override == TRUE){
                                            if("0" %in% unique(UMAP_df$cluster_ID)){
                                            print("Found unclustered boutons")
                                            color_length = length(unique(UMAP_df$cluster_ID))-1
                                            color_list <- scales::seq_gradient_pal("blue", "orange")(seq(0,1,length.out=color_length))
                                            print(color_list)
                                            color_list <- append(c('#828282'), color_list)
                                            } else {
                                            print("no boutons were unclustered")
                                            color_length = length(unique(UMAP_df$cluster_ID))
                                            color_list <- scales::seq_gradient_pal("blue", "orange")(seq(0,1,length.out=color_length))
                                            print(color_list)
                                            

                                            }

                                            #print(color_list)


                            }

                            umap_plot<-ggplot(UMAP_df, aes_string("x", "y",colour=colour_var)) +
                                        geom_point(alpha=0.8,size=3)+
                                        {if(scale_bool == FALSE)scale_colour_gradient(low='blue',high='red')}+
                                        {if(scale_bool)scale_colour_manual(values = color_list)}+
                                        #scale_size_continuous(name="Normalized amplitude, 4 mM Ca2+")+
                                        labs( x="UMAP 1",
                                              y="UMAP 2",
                                              title = plot_title,
                                             tag=tag_var)+
                                        coord_cartesian(xlim=c(-7,7), ylim=c(-7,7))+
                                        #coord_cartesian(xlim=c(x_min*1.2,x_max*1.2),ylim=c(y_min*1.2,y_max*1.2))+
                                        #guides(colour = guide_colourbar(position = "inside"))+
                                        scale_y_continuous(breaks=c(-5,0,5))+
                                        scale_x_continuous(breaks=c(-5,0,5))+
                                        guides(colour = guide_legend(reverse=FALSE, title = "Cluster ID"))+
                                        theme_tufte()+
                                        my.theme+
                                        theme(legend.position.inside=c(0.9,0.9),
                                              legend.justification.inside = "right",
                                              legend.title.align = 1,
                                              legend.margin = margin(0,0,0,0),

                                                legend.title = element_text(colour="black", size=34*scalefactor, family="sans",margin=margin(b=10)),
                                                legend.spacing.y = unit(2,'cm'))
                                                         

                                            if(global_output == TRUE){
                                                                nam <- plot_filename
                                                                assign(nam, umap_plot,envir = .GlobalEnv)
                                                                }
                                
                            save_plot(filename = paste0(plot_filename), plot=umap_plot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
color_list             

}

#this function simply adds the cluster annotation to the 






#indicate which columns do not need to be normalized between 0 and 1 (e.g. columns that already range from 0 to 1, like PiGlu)
exclude_cols<- c('PiGlu_05Ca','PiGlu_1Ca','PiGlu_2Ca','PiGlu_4Ca')
UMAP_corr_stats<- generate_UMAP(corr_stats,exclude_cols)



#indicate which columns do not need to be normalized between 0 and 1
exclude_cols<- c('PiGlu_05Ca','PiGlu_1Ca','PiGlu_2Ca','PiGlu_4Ca','Pr_05Ca','Pr_1Ca','Pr_2Ca','Pr_4Ca','RRP_ratio_05Ca','RRP_ratio_1Ca','RRP_ratio_2Ca','RRP_ratio_4Ca','N_sites')#, 'Q')
UMAP_corr_stats_w_MV<- generate_UMAP(corr_stats_w_MV, exclude_cols)

#cluster outputs with hdbscan
UMAP_corr_clustered= cluster_UMAP(UMAP_corr_stats)
UMAP_corr_w_MV_clustered= cluster_UMAP(UMAP_corr_stats_w_MV)

UMAP_color_list = plot_UMAP(UMAP_df = UMAP_corr_clustered, 
            colour_var = "cluster_ID", 
            plot_title = expression(atop("Clustered UMAP", "of bouton function") ), 
            plot_filename = "clustered_UMAP",
            tag_var = "E",
            global_output = TRUE,
            scale_bool = TRUE,
            colour_override = TRUE) 




UMAP_MV_color_list = plot_UMAP(UMAP_df = UMAP_corr_w_MV_clustered, 
            colour_var = "cluster_ID", 
            plot_title = str_wrap("Clustered UMAP with MV analysis",25), 
            plot_filename = "clustered_UMAP_w_MV", 
            tag_var = "E",
            global_output = FALSE,
            scale_bool = TRUE,
            colour_override = TRUE) 
print(UMAP_color_list)



#### strategy for Figure 4

# ABCD = original ABCD


#A = releaseProb_plot <- PiGlu
#B = CV_plot<- amplitude
#C = CV_plot<- tau_decay
#D = CV_plot<- dt

#E = UMAP clustered

#F = mean_quanta across Ca per clusterID
#G = 




# ###### 1) PLOT RELEASE PROBABILITY PER ROI  #######

                    releaseProb_calc$groupfact <- factor(releaseProb_calc$Ca_mM)    
                    releaseProb_calc$Ca<- as.factor(as.character( (releaseProb_calc$Ca) ) )
                    levels(releaseProb_calc$Ca) <- new_levels

                        

                    releaseProbPlot<-ggplot(releaseProb_calc, aes(x=Ca_mM, y=releaseProbability_perROI))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=1.2,alpha=0.35)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white", alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression(P[iGlu]),
                                    title="",
                                    tag="A"
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                            #guides(colour="none")+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
                            coord_cartesian(ylim=c(0,1.1),xlim=c(0.25,4.25))+
                            guides(colour = guide_legend(override.aes = list(colour = color_override,
                                                                                size = c(4,4,4,4),
                                                                                alpha=c(1,1,1,1)
                                                                                )
                                                                            )
                                                                        )+                    
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "right")


            save_plot(filename = "releaseProbPlot_v1", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                peak_probability<<- releaseProbPlot






# Function to create CV plot
create_CV_plot <- function(data, y_variable, y_label, tag_var, filename,y_upper_bound=1.75) {
  print(paste0("creating a CV violin plot for: ", y_variable) )

  x_upper_bound = 4 *1.5
  jitter_x_range = x_upper_bound/20
  if(y_upper_bound == 1.75){
    break_vec = c(0, 0.5, 1.0,1.5)
  } else {
    break_vec = c(25,50,75,100)
  }


  plot <- ggplot(data, aes(x = groupfact, y = get(y_variable), group = groupfact, colour = Ca)) +
    geom_sina(aes(group=groupfact, colour=groupfact),size=1.5,alpha=0.5)+
    geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
    

    labs(
      y = y_label,
      x = expression("["*Ca^'2+'*"] (mM)"),
      colour = "",
      title = "",
      tag = tag_var
    ) +
    scale_y_continuous(breaks =  c(0, 0.5, 1.0,1.5)) + #expand=c(0,0)
    coord_cartesian(ylim = c(0, y_upper_bound)) +
    {if (color_switch == FALSE) scale_colour_brewer(palette = "Dark2")} +
    {if (color_switch) scale_colour_manual(labels = new_levels, values = color_override)} +
    theme_tufte() +
    my.theme +
    theme(legend.position = "none")

  save_plot(filename = paste0(filename, ""), plot = plot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

    if(!is.null(tag_var)){assign(paste0(y_variable), plot, envir = .GlobalEnv)}
}



# # Create CV plots




create_CV_plot(releaseProb_calc, "CV_amplitude", expression("CV of Peak "*Delta*"F/F"), "B", "CV_amplitude_violin",y_upper_bound = 1.1)
create_CV_plot(releaseProb_calc, "CV_tau_decay", expression("CV of "*tau[decay]), "C", "CV_tau_violin",y_upper_bound=1.1)
create_CV_plot(releaseProb_calc, "CV_dt", expression("CV of "*Delta*"t"), "D", "CV_dt_violin",y_upper_bound=1.1)
#create_CV_plot(releaseProb_calc, "CV_t_half", expression("CV of "*t[1/2]), NULL, "CV_t_half_violin")
#create_CV_plot(releaseProb_calc, "CV_t_rise", expression("CV of "*t[rise]), NULL, "CV_t_rise_violin")
# # create_CV_plot(releaseProb_calc, "CV_t_decay", expression("CV of "*t[decay]), NULL, "CV_t_decay_violin")
# create_CV_plot(good_N_sites, "P_rw", expression("Release Probability, "*P[r]), NULL, "binomial_Pr")
# create_CV_plot(good_N_sites, "Q_w", expression("Quantal size, Q"), NULL, "binomial_Q")
# create_CV_plot(good_N_sites, "N_min", expression("No. Release Sites, "*N[min]), NULL, "binomial_N", 100)








###### Segregating data by UMAP clusters


#UMAP_clusterIDs<- UMAP_corr_clustered %>%  
UMAP_cluster_ID<- UMAP_corr_clustered$cluster_ID

corr_stats_clusterID<- corr_stats %>% bind_cols(clusterID = UMAP_cluster_ID) %>% select(trackROI_key, clusterID) 
clean_peaks_normTime_clustered<- left_join(clean_peaks_normTime, corr_stats_clusterID) %>% dplyr::filter(clusterID > 0 ) 


releaseProb_calc_clustered<-left_join(releaseProb_calc,corr_stats_clusterID) %>% dplyr::filter(clusterID != 0)
#coerce mean_quanta to 0 to get a better definition of mean_quanta)
releaseProb_calc_clustered$mean_quanta[is.na(releaseProb_calc_clustered$mean_quanta)] <- 0
releaseProb_calc_clustered$CV_amplitude[is.na(releaseProb_calc_clustered$CV_amplitude)] <- 0
releaseProb_calc_clustered$CV_tau_decay[is.na(releaseProb_calc_clustered$CV_tau_decay)] <- 0
releaseProb_calc_clustered$CV_dt[is.na(releaseProb_calc_clustered$CV_dt)] <- 0
releaseProb_calc_clustered$releaseProbability_perROI[is.na(releaseProb_calc_clustered$releaseProbability_perROI)] <- 0



MV_params_clustered<- left_join(MV_params, corr_stats_clusterID) %>% dplyr::filter(clusterID > 0) 
MV_params_clustered_4Ca<- MV_params_clustered %>% dplyr::filter( Ca_mM == 4)

clustered_MV_summary<- MV_params_clustered %>% ungroup() %>% 
                                               group_by(clusterID, Ca, Ca_mM) %>%
                                               summarise(n_obs = n(), 
                                                            avg_N_sites = mean(N),
                                                            sd_N_sites = sd(N),
                                                            se_N_sites = sd_N_sites/sqrt(n_obs),

                                                            mean_RRP_ratio = mean(RRP_ratio),
                                                            sd_RRP_ratio = sd(RRP_ratio),
                                                            se_RRP_ratio = sd_RRP_ratio/sqrt(n_obs),
                                               
                                                            mean_Pr = mean(P_rw),
                                                            sd_Pr = sd(P_rw),
                                                            se_Pr = sd_Pr/sqrt(n_obs),

                                                            pop_mean_quanta = mean(mean_quanta),
                                                            sd_quanta = sd(mean_quanta),
                                                            se_quanta = sd_quanta/sqrt(n_obs),


                                                           mean_Qw = mean(Q_w),
                                                           sd_Qw = sd(Q_w),
                                                           se_Qw = sd_Qw/sqrt(n_obs),
                                                            )

   
clustered_releaseProb_summary <- releaseProb_calc_clustered %>% group_by(clusterID, Ca,Ca_mM) %>% 
                                                            summarise(n_obs = n(),
                                                                    mean_P_iGlu = mean(releaseProbability_perROI,na.rm=TRUE),
                                                                    se_P_iGlu = sd(releaseProbability_perROI,na.rm=TRUE)/sqrt(n_obs),

                                                                    mean_mean_quanta = mean(mean_quanta,na.rm=TRUE),
                                                                    se_mean_quanta = sd(mean_quanta,na.rm=TRUE)/sqrt(n_obs),

                                                                    mean_CV_amplitude = mean(CV_amplitude,na.rm=TRUE),
                                                                    se_CV_amplitude = sd(CV_amplitude,na.rm=TRUE)/sqrt(n_obs),
                                                                    
                                                                    mean_CV_tau_decay = mean(CV_tau_decay, na.rm=TRUE),
                                                                    se_CV_tau_decay = sd(CV_tau_decay,na.rm=TRUE)/sqrt(n_obs),
                                                                    
                                                                    mean_CV_dt = mean(CV_dt, na.rm=TRUE),
                                                                    se_CV_dt = sd(CV_dt,na.rm=TRUE)/sqrt(n_obs),
                                                                    
                                                                    
                                                                    mean_mean_dt = mean(mean_dt,na.rm=TRUE),
                                                                    se_dt = sd(mean_dt,na.rm=TRUE)/sqrt(n_obs),
                                                                    
                                                                    

                                                                    mean_mean_tau_decay = mean(mean_tau_decay,na.rm=TRUE),
                                                                    se_tau_decay = sd(mean_tau_decay,na.rm=TRUE)/sqrt(n_obs),
                                                                    
                                                                    mean_mean_t_rise = mean(mean_t_rise,na.rm=TRUE),
                                                                    se_t_rise = sd(mean_t_rise,na.rm=TRUE)/sqrt(n_obs),

                                                                    
                                                                    mean_mean_t_decay = mean(mean_t_decay,na.rm=TRUE),
                                                                    se_t_decay = sd(mean_t_decay,na.rm=TRUE)/sqrt(n_obs),

                                                                    
                                                                    mean_mean_t_half = mean(mean_t_half,na.rm=TRUE),
                                                                    se_t_half = sd(mean_t_half,na.rm=TRUE)/sqrt(n_obs)
                                                                    )






SI_scatterplot_omit05Ca<-function(df, y_var, y_expression, filename, tag_var){

                            releaseProbPlot<-ggplot(df, aes_string(x="Ca_mM", y=y_var, colour = "clusterID"))+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=y_expression,
                                    #title=expression(atop("Bouton activity", 
                                    #                        "by cluster ID") ),
                                    tag=tag_var
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(limits=c(0.75,4.25),breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
                            coord_cartesian(xlim=c(1,4.25))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")
                            
                            save_plot(filename = filename, plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                            nam<- filename
                            assign(nam, releaseProbPlot,envir = .GlobalEnv)
}

SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "CV_amplitude", y_expression = expression("CV of Peak "*Delta*"F/F"), filename = "clustered_CV_amplitude", tag_var = "B" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "CV_tau_decay", y_expression = expression("CV of "*tau[decay]), filename = "clustered_CV_tau_decay", tag_var = "C" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "CV_dt", y_expression = expression("CV of "*Delta*"t"), filename = "clustered_CV_dt", tag_var = "D" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "mean_tau_decay", y_expression = expression(tau[decay]~"(ms)"), filename = "clustered_tau_decay", tag_var = "E" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "mean_dt", y_expression = expression(Delta*"t (ms)"), filename = "clustered_dt", tag_var = "F" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "mean_t_half", y_expression = expression(t[1/2]~"(ms)"), filename = "clustered_t_half", tag_var = "G" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "mean_t_rise", y_expression = expression(t[rise]~"(ms)"), filename = "clustered_t_rise", tag_var = "H" )
SI_scatterplot_omit05Ca(df = releaseProb_calc_clustered, y_var = "mean_t_decay", y_expression = expression(t[decay]~"(ms)"), filename = "clustered_t_decay", tag_var = "I" )








                    

                    releaseProbPlot<-ggplot(releaseProb_calc_clustered, aes(x=Ca_mM, y=mean_quanta, colour = clusterID))+
                            geom_hline(yintercept = 1, alpha=0.8,size=0.7,lty="dashed",colour="black")+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+ 
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression("E"[Delta*"F/F"] / "S"[Delta*"F/F"]),
                                    title=expression(atop("Estimated SV release", 
                                                            "by cluster ID") ),
                                    tag="G"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            coord_cartesian(ylim=c(0,11),xlim=c(0.25,4.25))+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(1,5,10))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "mean_quanta_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_quanta<<- releaseProbPlot


       

                    releaseProbPlot<-ggplot(releaseProb_calc_clustered, aes(x=Ca_mM, y=releaseProbability_perROI, colour = clusterID))+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression(P[iGlu]),
                                    title=expression(atop("Bouton activity", 
                                                            "by cluster ID") ),
                                    tag="F"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(0,0.5,1.0))+
                            coord_cartesian(ylim=c(0,1.1),xlim=c(0.25,4.25))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "PiGlu_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_PiGlu<<- releaseProbPlot

                





       releaseProbPlot<-ggplot(MV_params_clustered, aes(x=Ca_mM, y=P_rw, colour = clusterID))+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression("Uniform"~P[r]),
                                    title=expression('Release probability'),
                                    tag="J"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "Pr_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_Pr<<- releaseProbPlot





       releaseProbPlot<-ggplot(MV_params_clustered, aes(x=Ca_mM, y=RRP_ratio, colour = clusterID))+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression( "("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"/N[sites]),
                                    title=expression(atop("Fraction of RRP", "released per AP")),
                                    tag="K"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "RRP_ratio_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_RRP_ratio<<- releaseProbPlot



       releaseProbPlot<-ggplot(MV_params_clustered_4Ca, aes(x=as.numeric(clusterID), y=Q_w, colour = clusterID, group = clusterID))+
                            geom_sina(size=2,alpha=0.8)+
                            geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
    
                           
                            labs( x=expression("Cluster ID"),
                                    y=expression(Q),
                                    title=expression("Quantal size"),
                                    tag="L"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(breaks=c(1,2,3),labels=c("1","2","3"))+
                            scale_y_continuous(breaks=c(0.5,1.0,1.5))+
                            coord_cartesian(ylim = c(0.5,1.5))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "Q_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_Q<<- releaseProbPlot

print(unique(MV_params_clustered$clusterID))


       releaseProbPlot<-ggplot(MV_params_clustered_4Ca, aes(x=as.numeric(clusterID), y=N, colour = clusterID, group = clusterID))+
                            geom_sina(size=2,alpha=0.8)+
                            geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
    
                           
                            labs( x=expression("Cluster ID"),
                                    y=expression(N[sites]),
                                    title=expression(atop("RRP size", 
                                                            "by cluster ID") ),
                                    tag="H"
                                    )+
                            scale_colour_manual(values=UMAP_color_list)+
                            scale_x_continuous(breaks=c(1,2,3),labels=c("1","2","3"))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")


            save_plot(filename = "N_sites_clustered", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                clustered_N_sites<<- releaseProbPlot



    avgtrace_Plot<-ggplot(clean_peaks_normTime_clustered, aes(x=new_normTime, y = dFF,  group=interaction(Ca,clusterID), colour=clusterID))+
                                                    stat_summary_bin(aes(group=clusterID), geom="line", fun=mean, binwidth=interFrame.var,size=1.5,alpha=0.8)+
                                                    
                                                    geom_vline(xintercept = 0, lty="dashed",colour='grey21',size=0.9)+
                                                    
                                                    labs( x="Normalized time (s)",
                                                          y=expression(Delta*"F/F"),
                                                          title=expression("Averaged iGluSnFR3 response by cluster ID"),
                                                          #colour="Segmented by:",
                                                                        tag = "A")+
                                                    scale_colour_manual(values=UMAP_color_list)+
                                                                #scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
                                                    scale_x_continuous(labels = c("0", "0.5"),breaks = c(0,0.5),limits=c(-0.055,0.65))+
                                                    theme_tufte()+
                                                    #guides(colour="none",fill="none")+
                                                    coord_cartesian(ylim=c(-0.01,5),xlim=c(-0.055,0.55))+
                                                    my.theme+
                                                    facet_grid(~Ca, labeller = label_parsed)+
                                                    theme(panel.spacing = unit(0.25, "cm"),
                                                        legend.position = "right")

                                                    save_plot(filename = "avg_traces_clustered", plot=avgtrace_Plot, plot_width=plot_height*2.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
            
                                                     clustered_avg_traces<<- avgtrace_Plot





# # # #### GET ALL THE SPATIAL DATA ##### 



###### for loop for plotting spatial data
        png_dir = "Y:\\Sam/paper1_datasets/singleAP_v2/imageFiles/zstack/zstackOutputs_sliced_SQ"
        prefix_str = "^MAX_GluSnFR3_SyPhy_"
        suffix_str = "_0pt5Ca_zstack_lowpwr__FLATTENED.png"
        all_lut_str = "_0pt5Ca_zstack_lowpwr__FLATTENED_COLOR.png"
        dir_regex = "GluSnFR3(.*?)region\\d\\d_"


                
        print("Taking a crack at generating png plots")
        subset_df<- releaseProb_calc_clustered %>% dplyr::filter(vid_key %in% vid_keys_to_save)

        ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d") 
        vid_key_vec <- unique(subset_df$vid_key)
        vars_to_plot<- c("clusterID","mean_quanta")#c("CV_amplitude",'mean_quanta','releaseProbability_perROI')#, #'release_ratio')#'releaseProbability_perROI')#c("releaseProbability_perROI")
        save_bool=TRUE
        guide_bool = TRUE

        guide_limits = list(c(1,3),c(1,20))#list(c(0,0.5),c(5,20),c(0,1))#,c(0.8,2.0))
        guide_titles = c("Cluster ID","SVs released")#c(expression(atop(CV[Delta*"F/F"], "1"~Ca^'2+')),expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4"~Ca^'2+')),expression(atop(P[iGlu], "0.5"~Ca^'2+')))#, expression(frac(E/S["4Ca"], E/S["2Ca"])))    #expression(P[iGlu], "0.5 mM"~Ca^'2+'))# # expression(CV[tau*"decay"])#
        
        #expression("Avg."~E[Delta*"F/F"] / S[Delta*"F/F"])
        tmp_tag_var_override = c("A","B")#c("L","L","L")#c("E","F", "G")


        ###scale bar params! set the left-corner 
        user_x = 20
        user_y = 24.86




        for (k in 1:length(vid_key_vec)){
            
                vid_df<- releaseProb_calc_clustered %>% dplyr::filter(vid_key == vid_key_vec[k]) %>% mutate(fix_Ca = paste0(Ca_mM, "Ca") ) %>% dplyr::filter(clusterID > 0)
                Ca_vec <- unique(vid_df$fix_Ca)
                print(Ca_vec)
        
                max_gradient_val = max(vid_df$releaseProbability_perROI)

                
                  count =0
                                    tmp_plot = NULL
                                    Ca_df_1<- vid_df %>% dplyr::filter(fix_Ca == "1Ca")
                                    Ca_df_4<- vid_df %>% dplyr::filter(fix_Ca == "4Ca")
                                    Ca_df_05<- vid_df %>% dplyr::filter(fix_Ca == "0.5Ca")
                                    Ca_df_list = list(Ca_df_1,Ca_df_4,Ca_df_05)#, Ca_df_saturated)

                                    for(j in 1:length(vars_to_plot)){

                                        if(save_bool==TRUE) { 
                                            guide_bool = TRUE
                                            var_to_plot = vars_to_plot[j]
                                            tmp_guide_title = guide_titles[j]
                                            tmp_limits = guide_limits[[j]]
                                            tmp_tag = tmp_tag_var_override[j]
                                            if(j %in% c(1,2,3)){tmp_scalebar = TRUE} else {tmp_scalebar = FALSE} ####### CHANGED SCALEBAR RULE
                                            Ca_df = Ca_df_list[[j]]
                                            print(Ca_df)

                                           tmp_plot = png_plotter(df=Ca_df,
                                                                    png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                                                    vars_to_plot=var_to_plot, 
                                                                    save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
                                                                    ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar)
                                    }

                                    }

                                                                    

                
        }
                                    
                                        



# # ##### step 1. establish subset of the dataframe that can be used to plot traces and physiology



# ##### N_sites probably breaks here, should be careful to feed data in from MV_analysis after retooling


tmp_trace_data <- clean_df %>% dplyr::filter(trackROI_key %in% ROIs_to_save)#, !ROI_key %in% unique(double_check$ROI_key))
print(unique(tmp_trace_data))
tmp_timezero_df<- suppressMessages(tmp_trace_data %>%  group_by_at(groupers) %>% summarise(time_zero = unique(firstStim[!is.na(firstStim)] ) ) )
tmp_trace_data<- suppressMessages(left_join(tmp_trace_data,tmp_timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))

tracked_ROIs <- unique(tmp_trace_data$trackROI_key)
tmp_tag_vars <- list(c("C","D"),#c("J"),#c("M","G"),#c("H","K"), #trace
                    c("E","G"),#c("K"),#,"H"),#c("I", "L"), #releaseprob
                    c("F","H"),#c("L"),
                    c(""))#c("M"))#,"I"))#c("J","M"))#, #CV amplitude
                    #c("K", "N")) #amplitude
for (i in 1:length(tracked_ROIs)) {
                    count = i
                    print(paste0("loop iter on count: ", i))
                    #print("Generating data visualizations for trackROI_key:", ROI)
                    ROI_title<-gsub("GluSnFR3-dish\\d\\d-plate\\d\\d-region\\d\\d-","",tracked_ROIs[i])


                    tmp_original_stats<- peak_stat_df %>% 
                                         dplyr::filter(trackROI_key == tracked_ROIs[i])
                    tmp_summary_data <- releaseProb_calc %>%
                                        dplyr::filter(trackROI_key == tracked_ROIs[i])
                    tmp_ROIs_with_peaks<- ROIs_with_peaks %>% 
                                        dplyr::filter(trackROI_key == tracked_ROIs[i],
                                                        peak_positive > 0)
                    positive_ROI_list<- unique(tmp_ROIs_with_peaks$ROI_key)


                   # print(tmp_summary_data)

                    tmp_trace_subset<- tmp_trace_data %>% 
                                        dplyr::filter(trackROI_key == tracked_ROIs[i])#,
                                                        #ROI_key %in% positive_ROI_list)
                    tmp_trace_subset$Ca<- as.factor(as.character( (tmp_trace_subset$Ca) ) )
                    levels(tmp_trace_subset$Ca) <- new_levels


                        
                    derived_amplitudes<- left_join(tmp_trace_subset,light_spont_avg) %>% 
                                         group_by(Ca, Ca_mM, trackROI_key,ROI_key) %>%
                                         dplyr::filter(new_normTime > 0, new_normTime < 0.5, peakID != "NotPeak") %>%
                                         summarise(amplitude = max(dFF, na.rm=TRUE),
                                                    quanta_calc = amplitude/spont_amplitude) %>%
                                         dplyr::filter(trackROI_key == tracked_ROIs[i])
                    derived_amplitudes$groupfact <- factor(derived_amplitudes$Ca_mM)    

                    max_amp = derived_amplitudes %>% dplyr::filter(Ca_mM == 4)
                    max_amp = max(max_amp$amplitude)
                    
                    print(paste0("The maximum amplitude for :",unique(derived_amplitudes$trackROI_key), "was ", max_amp," DeltaF/F"))


                    tmp_amplitude<-ggplot(derived_amplitudes, aes(x=Ca_mM, y=quanta_calc))+#linetype=segmentation,
                                             #geom_sina(aes(group=groupfact, colour=groupfact),size=2,alpha=0.5)+
                            
                                             geom_point(aes(group=groupfact,colour=groupfact),alpha=0.5, size=2)+
                                             
                                             geom_hline(yintercept = 1, size=0.7,lty="dashed",colour="black")+
                                             stat_summary(geom="line", fun.y = mean, size=1, alpha=0.85,colour="black" )+
                                             
                                             stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=0.8, alpha=1, colour="black")+
                                             stat_summary(geom="point", fun.y=mean, color="black",shape=21,stroke=1.1, size=2.5,fill="white",alpha=1)+
                                             
                                             labs( y=expression(E[Delta*"F/F"] / S[Delta*"F/F"]),
                                                    x=expression("["*Ca^'2+'*"] (mM)"),
                                                    title=expression(atop("Estimated SVs","released per AP")),
                                                    colour="",
                                                    tag = tmp_tag_vars[[3]][i]
                                                    )+
                                            {if(color_switch)scale_colour_manual(values=c("0.5" = color_override[1],"1" = color_override[2],"2" = color_override[3],"4" = color_override[4]) )}+
                                            #{if(color_switch)scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))}+
                                                    
                                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                                            scale_y_continuous(breaks=c(1,10,20,30))+
                                            coord_cartesian(ylim=c(0,30),xlim=c(0.5,4))+
                            
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position = "none",
                                                  plot.title=element_text(size=scalefactor*34, family="sans",hjust=0.5) 
                 
                                                  )
                    
                    nam <- paste("tmp_amplitude", count, sep = "")
                    assign(nam, tmp_amplitude,envir = .GlobalEnv)
                    
                      

                    # tmp_CV_amplitude<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_amplitude))+#linetype=segmentation,
                    #                           geom_line(size=1,color="black",alpha=0.85)+
                    #                             geom_point(color="black",shape=21,stroke=1.1, size=2.5,fill="white")+
                                                                 
                    #                             labs( x=expression("["*Ca^'2+'*"] (mM)"),
                    #                                     y=expression("CV"[Delta*F/F]),
                    #                                     title="",
                    #                                     tag = tmp_tag_vars[[3]][i]
                    #                                     )+
                    #                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                    #                             {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                    #                             #guides(colour="none")+
                    #                             scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                    #                             scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
                    #                             coord_cartesian(ylim=c(0,0.5),xlim=c(0.25,4.25))+
                    #                                            #, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                    #                         theme_tufte()+
                    #                         my.theme+
                    #                         theme(legend.position = "none")

                    #     nam <- paste("tmp_CV_amplitude", count, sep = "")
                    #     assign(nam, tmp_CV_amplitude,envir = .GlobalEnv)
                                
                        


                    #     tmp_CV_tau_decay<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_tau_decay))+#linetype=segmentation,
                    #                          geom_line(size=2,color="black")+
                    #                          geom_point(size=4,color="black")+
                    #                          labs( y=expression("CV of "*tau[decay]*", ms"),
                    #                                 x=expression("[Ca"^{"2+"}*"], mM"),
                    #                                 title="",
                    #                                 colour=""
                    #                                 )+
                    #                         coord_cartesian(ylim=c(0,1.0))+
                    #                         scale_x_continuous(breaks=c(0.5,1,2,4))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                    #                         theme_tufte()+
                    #                         my.theme+
                    #                         theme(legend.position = "none")
                                            
                       



                    #     tmp_CV_dt<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_dt))+#linetype=segmentation,
                    #                          geom_line(size=2,color="black")+
                    #                          geom_point(size=4,color="black")+
                    #                          labs( y=expression("CV of "*Delta*"t, ms"),
                    #                                 x=expression("[Ca"^{"2+"}*"], mM"),
                    #                                 title="",
                    #                                 colour=""
                    #                                 )+
                    #                         coord_cartesian(ylim=c(0,1.0))+
                    #                         scale_x_continuous(breaks=c(0.5,1,2,4))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                    #                         theme_tufte()+
                    #                         my.theme+
                    #                         theme(legend.position = "none")
                                            
                      


                            tmp_fit_predict<- MV_analysis$fit_predict %>% dplyr::filter(trackROI_key == tracked_ROIs[i])
                            tmp_N_sites <- MV_analysis$clean_N_sites %>% dplyr::filter(trackROI_key == tracked_ROIs[i]) 

                          
                      

                            if(max(unique(tmp_N_sites$mean_quanta), na.rm=TRUE) < 15){
                            x_limit = 20
                                
                            } else {
                                x_limit = 30
                            }
                            y_limit = 1.5*max(unique(tmp_fit_predict$var_quanta),na.rm=TRUE)

                            #P1_nonunif = tmp_nonunif_fit_values %>% dplyr::filter(Ca == "1Ca")
                            P05 = tmp_N_sites %>% dplyr::filter(Ca == "0pt5Ca")
                            P1 = tmp_N_sites %>% dplyr::filter(Ca == "1Ca")
                            P2 = tmp_N_sites %>% dplyr::filter(Ca == "2Ca")
                            P4 = tmp_N_sites %>% dplyr::filter(Ca == "4Ca")

                                        tmp_N_sites_plot<-ggplot(tmp_fit_predict, aes(x=mean_quanta, y=var_quanta, group=trackROI_key))+ 
                                                                                    
                                                                    geom_line( colour="black", size=0.9, alpha=0.9, lty="solid")+
                                                                    
                                                                    geom_point(data=tmp_N_sites,aes(x=mean_quanta,y=var_quanta, colour=Ca), size=6,alpha=1)+
                                                                    
                                                                    
                                                                    #unif
                                                                    geom_text(label = paste0("N_sites = ", round(unique(tmp_N_sites$N),2 )), size=7, x= x_limit*0.5, y=y_limit*0.85) +
                                                                    
                                                                    labs(x=expression(bar(x)~"("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"),
                                                                            y=expression(sigma^2~"("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"),
                                                                            tag=tmp_tag_vars[[4]][i],
                                                                            title = expression(atop("Mean-variance", "analysis"))
                                                                            #title="Mean-Variance plot with binomial fit"#,
                                                                            #subtitle=current_ROI
                                                                            )+
                                                                    scale_colour_manual(values = color_override)+
                                                                    scale_x_continuous(breaks=c(0,5,10,15,20,30))+
                                                                    scale_y_continuous()+     
                                                                    coord_cartesian(xlim=c(0,x_limit),ylim=c(0,y_limit))+
                                                                    
                                                                    #facet_grid(~trackROI_key)+
                                                                    
                                                                    theme_tufte()+
                                                                    my.theme+
                                                                    theme(legend.position = "none",
                                                                            plot.margin = unit(c(0,0,0,0),"cm"),
                                                                            plot.title=element_text(size=scalefactor*34, family="sans",hjust=0.5))

                                                                    save_plot(filename = paste0("mean-variance_nlsFits_", tracked_ROIs[i]), plot=tmp_N_sites_plot, plot_width=plot_height*2,plot_height=plot_height*1.2, scale_f = scalefactor,dpi_val = 300)


                       
                            
                            nam <- paste("tmp_N_sites", count, sep = "")
                            assign(nam, tmp_N_sites_plot,envir = .GlobalEnv)
                            
                         #  


                        tmp_releaseProb<- ggplot(tmp_summary_data, aes(x=Ca_mM, y=releaseProbability_perROI, colour="Observed PiGlu"))+
                                                
                                                geom_line(size=1,alpha=0.85)+
                                                geom_point(shape=21,stroke=1.1, size=2.5,fill="white")+
                                                geom_line(data= tmp_N_sites,aes(x=Ca_mM, y=P_rw, colour="Uniform Pr"),size=1,alpha=0.85)+
                                                geom_point(data= tmp_N_sites,aes(x=Ca_mM, y=P_rw, colour="Uniform Pr"),shape=21,stroke=1.1, size=2.5,fill="white")+                 
                                                labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                                        y=expression(P),
                                                        #title="",
                                                        tag = tmp_tag_vars[[2]][i]
                                                        )+
                                                scale_colour_manual(name = NULL, values=c("black","blue"), breaks=c("Observed PiGlu", "Uniform Pr"))+
                                                #{if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                                #{if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                guides(colour=guide_legend(nrow=2,byrow=TRUE))+
                                                scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                                                scale_y_continuous(breaks=c(0,0.5,1.0))+
                                                coord_cartesian(ylim=c(0,1.0),xlim=c(0.5,4))+
                                                                    
                                                theme_tufte()+
                                                my.theme+
                                                theme(legend.position = "top"#,
                                                        #legend.justification.inside =1  
                                                        #axis.text.x=element_text(colour="black", size=28*scalefactor, family="sans",angle=45, hjust=1)
                                                        )
                        nam <- paste("tmp_releaseProb", count, sep = "")
                        assign(nam, tmp_releaseProb,envir = .GlobalEnv)
                    


                        tmp_display_traces<- tmp_trace_subset #%>% dplyr::filter(!ROI_key %in% unique(double_check$ROI_key))


                        tmp_tracePlot<- ggplot(tmp_display_traces, aes(x=new_normTime, y = dFF,  group=1, colour=Ca,alpha=Ca))+
                                                    geom_vline(xintercept = 0, lty="longdash",colour='grey21',size=0.9)+
                                                    
                                                    geom_path(size=0.5)+
                                                    stat_summary_bin(aes(group=Ca), colour="black",geom="line", fun=mean, binwidth=interFrame.var,size=1,alpha=0.95)+
                                                    
                                                    labs( x="Normalized time (s)",
                                                          y=expression(Delta*"F/F"),
                                                          title=paste0("Bouton",gsub("ROI00", "  ", ROI_title)),
                                                          #colour=expression("[Ca"^{"2+"}*"], mM"),
                                                                        tag = tmp_tag_vars[[1]][i])+
                                                    {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                                    {if(color_switch)scale_colour_manual(labels=new_levels,values=color_override)}+
                                                    {if(color_switch)scale_alpha_manual(labels=new_levels,values=c(0.8,0.7,0.5,0.35))}+
                                                    coord_cartesian(ylim=c(0,12.5))+
                                                    scale_y_continuous(breaks=c(0,2.5,5,7.5,10,12.5))+
                                                    scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
                                                    theme_tufte()+
                                                    my.theme+
                                                    facet_grid(~Ca, labeller = label_parsed)+
                                                    theme(legend.title = element_text(colour="black", size=28, family="sans"),
                                                        legend.position="none",
                                                        strip.text.x = element_text(colour="black",size=30*scalefactor,family="sans"))

                                                        

                                                     ggsave(filename=paste0("avg_traces_", tracked_ROIs[i],"_"),plot=tmp_tracePlot, device="jpeg",dpi=600, bg="white",units="in",width=plot_height*1.8,height=plot_height)
                                                         
                                                       if(count <= 2){
                                                        nam <- paste("tmp_tracePlot", count, sep = "")
                                                        assign(nam, tmp_tracePlot,envir = .GlobalEnv)
                                                        }

                        # inset_traces = tmp_display_traces %>% dplyr::filter(Ca_mM == 0.5)
#                         tmp_inset<- ggplot(inset_traces, aes(x=new_normTime, y = dFF,  group=1, colour=Ca,alpha=Ca))+
#                                                     geom_vline(xintercept = 0, lty="longdash",colour='grey21',size=0.9)+
                                                    
#                                                     geom_path(size=0.5)+
#                                                     stat_summary_bin(aes(group=Ca), colour="black",geom="line", fun=mean, binwidth=interFrame.var,size=1,alpha=0.95)+
                                                    
#                                                     labs( x="Normalized time (s)",
#                                                           y=expression(Delta*"F/F")#,
#                                                           )+
#                                                     {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                                     {if(color_switch)scale_colour_manual(labels=new_levels,values=color_override)}+
#                                                     {if(color_switch)scale_alpha_manual(labels=new_levels,values=c(0.8,0.7,0.5,0.35))}+
#                                                     coord_cartesian(ylim=c(-0.1,2))+
#                                                     scale_y_continuous(breaks=c(0,1,2))+
#                                                     scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
#                                                     theme_tufte()+
#                                                     my.theme+
#                                                     #facet_grid(~Ca, labeller = label_parsed)+
#                                                     theme(legend.title = element_text(colour="black", size=28, family="sans"),
#                                                         legend.position="none",
#                                                         strip.text.x = element_text(colour="black",size=34*scalefactor,family="sans"))

                                                        

# #                                                     ggsave(filename=paste0("inset_traces_05Ca_", tracked_ROIs[i],"_.jpeg"),plot=tmp_inset, device="jpeg",dpi=600, bg="white",units="in",width=plot_height*1.8,height=plot_height*)
#                                                     save_plot(filename=paste0("inset_traces_05Ca",tracked_ROIs[i],"_"),plot=tmp_inset, plot_width=28*0.5,plot_height=12, scale_f = scalefactor,dpi_val = 600)
                      


                        # gs = list(tmp_tracePlot,  #1
                        #          tmp_releaseProb,
                        #          tmp_N_sites,       #2
                        #          tmp_amplitude)  #3
                        # margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
                    



                        # hlay <- rbind(c(1,1,1,1,1,1,2,2,2,3,3,3),
                        #                c(1,1,1,1,1,1,2,2,2,3,3,3),
                        #                c(1,1,1,1,1,1,2,2,2,3,3,3),
                                       
                        #                c(1,1,1,1,1,1,4,4,4,4,4,4),
                        #                c(1,1,1,1,1,1,4,4,4,4,4,4),
                        #                c(1,1,1,1,1,1,4,4,4,4,4,4)
                        #            )
                
                    

                    
                      #  ROI_report = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
                      #  save_plot(filename=paste0("ROI_report_",tracked_ROIs[i],"_"),plot=ROI_report, plot_width=28,plot_height=12, scale_f = scalefactor,dpi_val = 600)



}

                               

        





list.output = list(UMAP_corr = UMAP_corr_stats, UMAP_corr_MV = UMAP_corr_stats_w_MV, cluster_UMAP_corr = UMAP_corr_clustered, cluster_UMAP_corr_MV = UMAP_corr_w_MV_clustered,
                        MV_params = MV_params, MV_params_clustered = MV_params_clustered, releaseProb_calc_clustered = releaseProb_calc_clustered,
                        MV_params_clustered_4Ca = MV_params_clustered_4Ca,
                        clustered_MV_summary = clustered_MV_summary,
                        clustered_releaseProb_summary = clustered_releaseProb_summary)
list.output
}
						
