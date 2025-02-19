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



subDir <- "track_singleAP_phys_v6_umap"   


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



releaseProb_stats<- function(peak_stat_df = peak_stats, df = df_traces,  groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops,ROIs_to_save = ROIs_to_save){




                        #### establish variables that will be used for generating datasets


                        vid_keys_to_save = unique(str_extract(ROIs_to_save, "dish\\d\\d-plate\\d\\d-region\\d\\d") ) 
                        color_switch = !is.null(color_override)
                        plot_height = 8
                        interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE) 
                        


                    #    remove_neurons = c("dish05-plate03-region02", "dish05-plate06-region01", "dish05-plate06-region02")#,"dish05-plate01-region01")
                        #### clean up dataframe by merging with peakStats

                        peak_stat_df <- peak_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID), 
                                                                Ca_mM = case_when(Ca == "8Ca" ~ 8,
                                                                                  Ca == "0pt6Ca" ~ 0.6,
                                                                                  Ca == "0pt8Ca" ~ 0.8,
                                                                                  Ca == "1Ca" ~ 1,
                                                                                  Ca == "1pt2Ca" ~ 1.2,
                                                                                  Ca == "1pt6Ca" ~ 1.6,
                                                                                  Ca == "2Ca" ~ 2,
                                                                                  Ca == "4Ca" ~ 4),
                                                                imaging_region = gsub("-repl\\d\\d", "", vid_key)) #%>%
                                                 # dplyr::filter(!imaging_region %in% remove_neurons)
                        
                        light_df <- df[,!(names(df) %in% remove_col)] %>% mutate(Ca_mM = case_when(Ca == "8Ca" ~ 8,
                                                                                  Ca == "0pt6Ca" ~ 0.6,
                                                                                  Ca == "0pt8Ca" ~ 0.8,
                                                                                  Ca == "1Ca" ~ 1,
                                                                                  Ca == "1pt2Ca" ~ 1.2,
                                                                                  Ca == "1pt6Ca" ~ 1.6,
                                                                                  Ca == "2Ca" ~ 2,
                                                                                  Ca == "4Ca" ~ 4))

                       
                        ROIs_to_remove <- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are finite
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%
                                                 dplyr::filter( is.na(dFF))
                        ROIs_to_remove_timeskip<- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are do not have a frame skip in time of interest
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%            
                                                 summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
                                                 dplyr::filter(maxInterframe > 0.18) %>%
                                                 select(ROI_key,maxInterframe) 

#print(double_check)
                        #dplyr::filter(maxInterframe > 0.1)
print(paste0("These are the bad recordings: ", unique(ROIs_to_remove$ROI_key)))
#print(paste0("These are the bad recordings with a frame jump during post-stimulus: ", unique(ROIs_to_remove_timeskip$ROI_key)))

                       
                        # print("These should be the borked ROIs we found (drift correction failures?)")
                        # print(unique(ROIs_to_remove$ROI_key))
                        # check<- ROIs_to_remove
                        ROIs_to_remove <- unique(ROIs_to_remove$ROI_key)     
                        ROIs_to_remove_timeskip<-unique(ROIs_to_remove_timeskip$ROI_key)                                             

                        clean_df <- light_df %>%  dplyr::filter(!ROI_key %in% ROIs_to_remove, !ROI_key %in% ROIs_to_remove_timeskip) %>% 
                                                  mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) #%>%
                                                  #dplyr::filter(!imaging_region %in% remove_neurons)
                        
                        replicate_groupers= c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","vid_key","trackROI_key")
                        colname_toFind = "exposeNum"
                        df_mutated<- find_replicates(df = clean_df, groupers=replicate_groupers, colname_toFind=colname_toFind)
                        clean_df<-df_mutated
                        print("These were the replicates we found after cleaning")
                        print(unique(clean_df$replicates))
                        rm(df_mutated)

                        print(unique(clean_df$Ca))

                        clean_peaks<- suppressMessages(left_join(clean_df,peak_stat_df) %>% dplyr::filter(windowedPeakID != "NotPeak") ) 
                        timezero_df<- suppressMessages(clean_peaks %>%  group_by_at(groupers) %>% summarise(time_zero = unique(firstStim[!is.na(firstStim)] ) ) )
                        clean_peaks_normTime<- suppressMessages(left_join(clean_peaks,timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))
                        

                        
                        new_levels<- c(expression("0.6 mM "*Ca^'2+'),expression("0.8 mM "*Ca^'2+'),expression("1 mM "*Ca^'2+'),
                                         expression("1.2 mM "*Ca^'2+'), expression("1.6 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'),
                                         expression("4 mM "*Ca^'2+'), expression("8 mM "*Ca^'2+') )  
                        
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










###### CALCULATE RELEASE PROBABILITY AND OTHER PARAMS PER ROI #######
                    

                    

                    ROIs_with_peaks<- suppressMessages(clean_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5), !ROI_key %in% ROIs_to_remove) %>%  
                                                                group_by(Ca, Ca_mM, vid_key,trackROI_key,ROI_key,ROINumber,exposeNum,replicates) %>%  
                                                                summarise(peak_detected = length(unique(peakID))-1,
                                                                                peak_positive = case_when(peak_detected > 0 ~ 1,
                                                                                                          peak_detected == 0 ~ 0) ) 
                                                                )

                    silent_ROIs <- ROIs_with_peaks %>% dplyr::filter(peak_positive == 0)
                    silent_ROIs_list <- unique(silent_ROIs$ROI_key)     
                        
                    true_noise_amplitude<- clean_df %>% dplyr::filter(ROI_key %in% silent_ROIs_list) %>% 
                                                        dplyr::filter(stimEpoch == 1, absoluteTime < (firstStim +0.2)) %>%
                                                        group_by(Ca, Ca_mM, vid_key,trackROI_key,ROI_key,ROINumber,exposeNum,replicates) %>%  
                                                        summarise(amplitude = max(dFF, na.rm=TRUE)) %>%
                                                        mutate(quanta_calc = amplitude/0.4,
                                                                 Ca_expr = case_when(Ca_mM == 0.6 ~ "0pt6Ca",
                                                                                     Ca_mM == 8 ~ "8Ca",
                                                                                                        Ca_mM == 0.8 ~ "0pt8Ca",
                                                                                                        Ca_mM == 1 ~ "1Ca",
                                                                                                        Ca_mM == 1.2 ~ "1pt2Ca",
                                                                                                        Ca_mM == 1.6 ~ "1pt6Ca",
                                                                                                        Ca_mM ==  2 ~ "2Ca",
                                                                                                        Ca_mM == 4 ~ "4Ca"))        
                    print("This is the true noise amplitude df")
                    print(true_noise_amplitude)
                    





                    droppers = c("sensor","marker","segmentation","TTL_start","dish","plate","region","peakID","protocol",'exposeNum', "windowedPeakID","Ca")
                    # cleanup_peak_stat_data_point <- peak_stat_df %>% dplyr::filter(ROI_key == "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005", interSpike_ms < 100)
                    # cleanup_peak_stat_df <- peak_stat_df %>% dplyr::filter(ROI_key != "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005" )
                    # cleanup_merge <- bind_rows(cleanup_peak_stat_df, cleanup_peak_stat_data_point)
                    cleanup <- peak_stat_df[,!(names(peak_stat_df) %in% droppers)]

                    # #spont_avg_data for quantal size calc
                    # light_spont_avg <- spont_avgs %>% mutate(Ca_mM = case_when(Ca == "\"0.5 mM \" * Ca^\"2+\"" ~ 0.5,
                    #                                                                 Ca == "\"1 mM \" * Ca^\"2+\"" ~  1,
                    #                                                                 Ca == "\"2 mM \" * Ca^\"2+\"" ~ 2,
                    #                                                                 Ca == "\"4 mM \" * Ca^\"2+\"" ~ 4),
                    #                                         spont_amplitude = mean_amplitude ) %>%
                    #                                 ungroup() %>%
                    #                                 select(Ca_mM,spont_amplitude)
                    # #print("Here's a glimpse of the spont_avg_df")
                    #print(head(light_spont_avg))

                    light_peak_stat_df <- suppressMessages(cleanup %>% 
                                                                            mutate(quanta_calc = amplitude/0.4, 
                                                                                    Ca_expr = case_when(Ca_mM == 0.6 ~ "0pt6Ca",
                                                                                                        Ca_mM == 0.8 ~ "0pt8Ca",
                                                                                                        Ca_mM == 1.0 ~ "1Ca",
                                                                                                        Ca_mM == 1.2 ~ "1pt2Ca",
                                                                                                        Ca_mM == 1.6 ~ "1pt6Ca",
                                                                                                        Ca_mM ==  2.0 ~ "2Ca",
                                                                                                        Ca_mM == 4.0 ~ "4Ca",
                                                                                                        Ca_mM == 8.0 ~ "8Ca")) )

                    #light_peak_and_noise_stat_df <- suppressMessages(full_join(light_peak_stat_df,true_noise_amplitude))
                    print("Here's a glimpse of the quanta_calc")
                    print(unique(light_peak_stat_df$quanta_calc))
                    light_peak_stat_df$Ca_expr <- factor(light_peak_stat_df$Ca_expr, levels = c("0pt6Ca","0pt8Ca","1Ca","1pt2Ca","1pt6Ca","2Ca","4Ca","8Ca"))
                    unique(light_peak_stat_df$Ca_expr)
                    levels(light_peak_stat_df$Ca_expr) <- new_levels

                                                            
                    

                    print("These are the columns that exist in peak_stat_df")
                    #print(names(peak_stat_df))
                    print(names(light_peak_stat_df))

                    print(head(light_peak_stat_df))
                         
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

#Generate RRP vs. PiGlu at 0.5 mM
#### mean_quanta
get_RRP_06Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_08Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_1Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_12Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_16Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_2Ca <- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_4Ca <- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)
get_RRP_8Ca <-releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,mean_quanta)

names(get_RRP_06Ca)[names(get_RRP_06Ca) == 'mean_quanta'] <- 'mean_quanta_06Ca'
names(get_RRP_08Ca)[names(get_RRP_08Ca) == 'mean_quanta'] <- 'mean_quanta_08Ca'
names(get_RRP_1Ca)[names(get_RRP_1Ca) == 'mean_quanta'] <- 'mean_quanta_1Ca'
names(get_RRP_12Ca)[names(get_RRP_12Ca) == 'mean_quanta'] <- 'mean_quanta_12Ca'
names(get_RRP_16Ca)[names(get_RRP_16Ca) == 'mean_quanta'] <- 'mean_quanta_16Ca'
names(get_RRP_2Ca)[names(get_RRP_2Ca) == 'mean_quanta'] <- 'mean_quanta_2Ca'
names(get_RRP_4Ca)[names(get_RRP_4Ca) == 'mean_quanta'] <- 'mean_quanta_4Ca'
names(get_RRP_8Ca)[names(get_RRP_8Ca) == 'mean_quanta'] <- 'mean_quanta_8Ca'


###mean amplitude
get_dFF_06Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_08Ca <-releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_1Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_12Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_16Ca <-releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_2Ca <- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_4Ca <- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)
get_dFF_8Ca <- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,mean_amplitude)


names(get_dFF_06Ca)[names(get_dFF_06Ca) == 'mean_amplitude'] <- 'mean_amplitude_06Ca'
names(get_dFF_08Ca)[names(get_dFF_08Ca) == 'mean_amplitude'] <- 'mean_amplitude_08Ca'
names(get_dFF_1Ca)[names(get_dFF_1Ca) == 'mean_amplitude'] <- 'mean_amplitude_1Ca'
names(get_dFF_12Ca)[names(get_dFF_12Ca) == 'mean_amplitude'] <- 'mean_amplitude_12Ca'
names(get_dFF_16Ca)[names(get_dFF_16Ca) == 'mean_amplitude'] <- 'mean_amplitude_16Ca'
names(get_dFF_2Ca)[names(get_dFF_2Ca) == 'mean_amplitude'] <- 'mean_amplitude_2Ca'
names(get_dFF_4Ca)[names(get_dFF_4Ca) == 'mean_amplitude'] <- 'mean_amplitude_4Ca'
names(get_dFF_8Ca)[names(get_dFF_8Ca) == 'mean_amplitude'] <- 'mean_amplitude_8Ca'


####CV_amplitude

get_CV_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)

get_CV_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)

get_CV_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)
get_CV_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,CV_amplitude)


names(get_CV_06Ca)[names(get_CV_06Ca) == 'CV_amplitude'] <- 'CV_amplitude_06Ca'
names(get_CV_08Ca)[names(get_CV_08Ca) == 'CV_amplitude'] <- 'CV_amplitude_08Ca'
names(get_CV_1Ca)[names(get_CV_1Ca) == 'CV_amplitude'] <- 'CV_amplitude_1Ca'
names(get_CV_12Ca)[names(get_CV_12Ca) == 'CV_amplitude'] <- 'CV_amplitude_12Ca'
names(get_CV_16Ca)[names(get_CV_16Ca) == 'CV_amplitude'] <- 'CV_amplitude_16Ca'
names(get_CV_2Ca)[names(get_CV_2Ca) == 'CV_amplitude'] <- 'CV_amplitude_2Ca'
names(get_CV_4Ca)[names(get_CV_4Ca) == 'CV_amplitude'] <- 'CV_amplitude_4Ca'
names(get_CV_8Ca)[names(get_CV_8Ca) == 'CV_amplitude'] <- 'CV_amplitude_8Ca'


#### var_amplitude
get_var_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)
get_var_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,var_amplitude)

names(get_var_06Ca)[names(get_var_06Ca) == 'var_amplitude'] <- 'var_amplitude_06Ca'
names(get_var_08Ca)[names(get_var_08Ca) == 'var_amplitude'] <- 'var_amplitude_08Ca'
names(get_var_1Ca)[names(get_var_1Ca) == 'var_amplitude'] <- 'var_amplitude_1Ca'
names(get_var_12Ca)[names(get_var_12Ca) == 'var_amplitude'] <- 'var_amplitude_12Ca'
names(get_var_16Ca)[names(get_var_16Ca) == 'var_amplitude'] <- 'var_amplitude_16Ca'
names(get_var_2Ca)[names(get_var_2Ca) == 'var_amplitude'] <- 'var_amplitude_2Ca'
names(get_var_4Ca)[names(get_var_4Ca) == 'var_amplitude'] <- 'var_amplitude_4Ca'
names(get_var_8Ca)[names(get_var_8Ca) == 'var_amplitude'] <- 'var_amplitude_8Ca'


#### PiGlu
get_PiGlu_06<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_08<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_1<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_12<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_16<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_2<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_4<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)
get_PiGlu_8<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,releaseProbability_perROI)

names(get_PiGlu_06)[names(get_PiGlu_06) == 'releaseProbability_perROI'] <- 'PiGlu_06Ca'
names(get_PiGlu_08)[names(get_PiGlu_08) == 'releaseProbability_perROI'] <- 'PiGlu_08Ca'
names(get_PiGlu_1)[names(get_PiGlu_1)   == 'releaseProbability_perROI'] <- 'PiGlu_1Ca'
names(get_PiGlu_12)[names(get_PiGlu_12) == 'releaseProbability_perROI'] <- 'PiGlu_12Ca'
names(get_PiGlu_16)[names(get_PiGlu_16) == 'releaseProbability_perROI'] <- 'PiGlu_16Ca'
names(get_PiGlu_2)[names(get_PiGlu_2)   == 'releaseProbability_perROI'] <- 'PiGlu_2Ca'
names(get_PiGlu_4)[names(get_PiGlu_4)   == 'releaseProbability_perROI'] <- 'PiGlu_4Ca'
names(get_PiGlu_8)[names(get_PiGlu_8)   == 'releaseProbability_perROI'] <- 'PiGlu_8Ca'


#####mean_tau_decay
get_mean_tau_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)
get_mean_tau_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)

get_mean_tau_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,mean_tau_decay)


names(get_mean_tau_06Ca)[names(get_mean_tau_06Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_06Ca'
names(get_mean_tau_08Ca)[names(get_mean_tau_08Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_08Ca'
names(get_mean_tau_1Ca)[names(get_mean_tau_1Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_1Ca'
names(get_mean_tau_12Ca)[names(get_mean_tau_12Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_12Ca'
names(get_mean_tau_16Ca)[names(get_mean_tau_16Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_16Ca'
names(get_mean_tau_2Ca)[names(get_mean_tau_2Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_2Ca'
names(get_mean_tau_4Ca)[names(get_mean_tau_4Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_4Ca'
names(get_mean_tau_8Ca)[names(get_mean_tau_8Ca) == 'mean_tau_decay'] <- 'mean_tau_decay_8Ca'




#####CV_tau_decay
get_CV_tau_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)
get_CV_tau_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,CV_tau_decay)


names(get_CV_tau_06Ca)[names(get_CV_tau_06Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_06Ca'
names(get_CV_tau_08Ca)[names(get_CV_tau_08Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_08Ca'
names(get_CV_tau_1Ca)[names(get_CV_tau_1Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_1Ca'
names(get_CV_tau_12Ca)[names(get_CV_tau_12Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_12Ca'
names(get_CV_tau_16Ca)[names(get_CV_tau_16Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_16Ca'
names(get_CV_tau_2Ca)[names(get_CV_tau_2Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_2Ca'
names(get_CV_tau_4Ca)[names(get_CV_tau_4Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_4Ca'
names(get_CV_tau_8Ca)[names(get_CV_tau_8Ca) == 'CV_tau_decay'] <- 'CV_tau_decay_8Ca'



####mean_dt

get_mean_dt_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)
get_mean_dt_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,mean_dt)


names(get_mean_dt_06Ca)[names(get_mean_dt_06Ca) == 'mean_dt'] <- 'mean_dt_06Ca'
names(get_mean_dt_08Ca)[names(get_mean_dt_08Ca) == 'mean_dt'] <- 'mean_dt_08Ca'
names(get_mean_dt_1Ca)[names(get_mean_dt_1Ca) == 'mean_dt'] <- 'mean_dt_1Ca'
names(get_mean_dt_12Ca)[names(get_mean_dt_12Ca) == 'mean_dt'] <- 'mean_dt_12Ca'
names(get_mean_dt_16Ca)[names(get_mean_dt_16Ca) == 'mean_dt'] <- 'mean_dt_16Ca'
names(get_mean_dt_2Ca)[names(get_mean_dt_2Ca) == 'mean_dt'] <- 'mean_dt_2Ca'
names(get_mean_dt_4Ca)[names(get_mean_dt_4Ca) == 'mean_dt'] <- 'mean_dt_4Ca'
names(get_mean_dt_8Ca)[names(get_mean_dt_8Ca) == 'mean_dt'] <- 'mean_dt_8Ca'



#####CV_dt
get_CV_dt_06Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_08Ca<- releaseProb_calc %>% dplyr::filter(Ca == "0pt8Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_1Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_12Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt2Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_16Ca<- releaseProb_calc %>% dplyr::filter(Ca == "1pt6Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)

get_CV_dt_2Ca<- releaseProb_calc %>% dplyr::filter(Ca == "2Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_4Ca<- releaseProb_calc %>% dplyr::filter(Ca == "4Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)
get_CV_dt_8Ca<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca") %>% ungroup() %>% select(trackROI_key,CV_dt)


names(get_CV_dt_06Ca)[names(get_CV_dt_06Ca) == 'CV_dt'] <- 'CV_dt_06Ca'
names(get_CV_dt_08Ca)[names(get_CV_dt_08Ca) == 'CV_dt'] <- 'CV_dt_08Ca'
names(get_CV_dt_1Ca)[names(get_CV_dt_1Ca) == 'CV_dt'] <- 'CV_dt_1Ca'
names(get_CV_dt_12Ca)[names(get_CV_dt_12Ca) == 'CV_dt'] <- 'CV_dt_12Ca'
names(get_CV_dt_16Ca)[names(get_CV_dt_16Ca) == 'CV_dt'] <- 'CV_dt_16Ca'
names(get_CV_dt_2Ca)[names(get_CV_dt_2Ca) == 'CV_dt'] <- 'CV_dt_2Ca'
names(get_CV_dt_4Ca)[names(get_CV_dt_4Ca) == 'CV_dt'] <- 'CV_dt_4Ca'
names(get_CV_dt_8Ca)[names(get_CV_dt_8Ca) == 'CV_dt'] <- 'CV_dt_8Ca'









stat_list<- list(get_RRP_06Ca,get_RRP_08Ca,get_RRP_1Ca, get_RRP_12Ca,  
                 get_RRP_16Ca, get_RRP_2Ca, get_RRP_4Ca,get_RRP_8Ca,  
                 
                 get_dFF_06Ca,get_dFF_08Ca,get_dFF_1Ca, get_dFF_12Ca,  
                 get_dFF_16Ca, get_dFF_2Ca, get_dFF_4Ca,get_dFF_8Ca,  
                 #get_var_05Ca, 
                 get_var_06Ca,get_var_08Ca,get_var_1Ca, get_var_12Ca,  
                 get_var_16Ca, get_var_2Ca, get_var_4Ca,get_var_8Ca,  
                 
                
                get_mean_dt_06Ca, get_mean_dt_08Ca, get_mean_dt_1Ca, get_mean_dt_12Ca,
                get_mean_dt_16Ca,get_mean_dt_2Ca, get_mean_dt_4Ca, get_mean_dt_8Ca, 
              
                get_mean_tau_06Ca, get_mean_tau_08Ca, get_mean_tau_1Ca, get_mean_tau_12Ca,
                get_mean_tau_16Ca,get_mean_tau_2Ca, get_mean_tau_4Ca, get_mean_tau_8Ca, 
              
                get_PiGlu_06, get_PiGlu_08, get_PiGlu_1, get_PiGlu_12,
                get_PiGlu_16,get_PiGlu_2, get_PiGlu_4, get_PiGlu_8, 
                
                get_CV_06Ca,get_CV_08Ca,get_CV_1Ca, get_CV_12Ca,
                get_CV_16Ca,get_CV_2Ca, get_CV_4Ca,get_CV_8Ca, 
                 
                get_CV_dt_06Ca, get_CV_dt_08Ca, get_CV_dt_1Ca, get_CV_dt_12Ca, 
                get_CV_dt_16Ca, get_CV_dt_2Ca, get_CV_dt_4Ca, get_CV_dt_8Ca,

                get_CV_tau_06Ca, get_CV_tau_08Ca,get_CV_tau_1Ca, get_CV_tau_12Ca,
                get_CV_tau_16Ca, get_CV_tau_2Ca, get_CV_tau_4Ca, get_CV_tau_8Ca
                )

corr_stats<- stat_list %>% reduce(full_join, by='trackROI_key')


# corr_stats<- corr_stats %>% mutate(RRP_score = case_when( mean_quanta_4Ca >= quantile(mean_quanta_4Ca, 0.50, na.rm=TRUE) ~ "red", # ##mean_quanta_4Ca - mean_quanta_2Ca > 0
#                                                                                    mean_quanta_4Ca <= quantile(mean_quanta_4Ca, 0.50, na.rm=TRUE) ~ "blue"), #mean_quanta_4Ca/mean_quanta_2Ca ## #mean_quanta_4Ca - mean_quanta_2Ca <= 0
#                                    release_ratio = mean_quanta_4Ca / mean_quanta_2Ca
#                                                                                  ) %>%
#                         dplyr::filter(!is.na(RRP_score))

### UMAP code starting point
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
standardize <- function(x){(x-min(x))/(max(x)-min(x))}
data_for_umap <- corr_stats %>% select(-trackROI_key)

data_for_umap[is.na(data_for_umap)] <- 0
data_for_umap[is.nan(data_for_umap)] <- 0
print(data_for_umap)

colrange = ncol(data_for_umap)
standardized_umap<- data_for_umap %>% mutate_at(1:colrange, funs(standardize(.))) # as.data.frame(scale(data_for_umap))
print(standardized_umap)

set.seed(125)


synapse.umap <- umap(standardized_umap)
df_xy<- as.data.frame(synapse.umap$layout)
colnames(df_xy)[1]<- "x"
colnames(df_xy)[2]<- "y"
df_xy<- df_xy %>% mutate(index = row_number())
df_data<- as.data.frame(synapse.umap$data)
df_data<-df_data %>% mutate(index = row_number())
df_annotated<-join(df_xy,df_data)



umap_plot<-ggplot(df_annotated, aes(x, y,colour=mean_amplitude_06Ca,size=mean_amplitude_8Ca)) +
            geom_point(alpha=0.5)+
            scale_colour_gradient(low='blue',high='red')+
            labs( x="UMAP 1",
                  y="UMAP 2",
                 tag="X")+
            #coord_cartesian(xlim=c(-4,4),ylim=c(-5,8))+
            #guides(colour='none',size='none')+
            scale_y_continuous(breaks=c(-8,-4,0,4,8))+
            scale_x_continuous(breaks=c(-8,-4,0,4,8))+
            theme_tufte()+
            my.theme+
            theme(legend.title = element_text(colour="black", size=34*scalefactor, family="sans"))
                             
    
ggsave(filename=paste0("umap_plot_nolegend_V1.jpeg"),plot=umap_plot, device="jpeg",dpi=600, units="in",width=plot_height*1.5,height=plot_height*1)

write.csv(df_annotated, file = "normalized_output_w_UMAP.csv")
write.csv(standardized_umap, file = "normalized_output_wo_UMAP.csv")
write.csv(data_for_umap, file = "raw_output.csv")
write.csv(corr_stats, file = "raw_output_with_identifiers.csv")

# check.data<- output_data$corr_stats %>% select(-RRP_score, -trackROI_key)
# synapse.umap <- umap(check.data_amplitude)
# check.data_amplitude[is.na(check.data_amplitude)] <- 0
# check.data_amplitude[is.nan(check.data_amplitude)] <- 0
#is.nan.data.frame <- function(x)
#do.call(cbind, lapply(x, is.nan))

#data123[is.nan(data123)] <- 0


# ggplot(df_annotated, aes(x, y, colour=cluster_label,group=cluster_label)) +
#     geom_point(alpha=0.5)+
#     geom_point(data=df_annotated_means, aes(x=x_hat, y=y_hat), colour='black',size=3)


# merge_labels<- corr_stats %>% select(trackROI_key,RRP_score,release_ratio)

# vis_distributions <- left_join(releaseProb_calc, merge_labels) %>% dplyr::filter(!is.na(RRP_score))


#                     vis_distributions$groupfact <- factor(vis_distributions$Ca_mM)    
#                     #vis_distributions$Ca<- as.factor(as.character( (vis_distributions$Ca) ) )
#                     #levels(vis_distributions$Ca) <- new_levels


# create_distribution_plots<- function(data, x_variable, y_variable, x_label,  y_label, group_var,filename,tag_var,add_text_labels,global_output){ 


#                 if(y_variable == 'mean_quanta'){
#                     add_hline = TRUE
#                 } else {
#                     add_hline = FALSE
#                 }
#                 releaseProbPlot<-ggplot(data, aes_string(x=x_variable, y=y_variable, group=group_var, colour=group_var))+
#                             #geom_sina(aes(group=groupfact),size=1.2,alpha=0.35)+
#                             stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85)+
                            
#                             {if(add_hline)geom_hline(yintercept = 1, size=0.7,lty="dashed",colour="black")}+
#                             #geom_smooth(formula=y~x, method='loess',alpha=1,se=FALSE,size=3,colour="black")+              
#                             stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1)+
#                             stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white",,alpha=1)+
                            
#                             labs( x=x_label,
#                                     y=y_label,
#                                     tag=tag_var
#                                     )+
#                             #guides(colour="none")+
#                             scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
#                             #scale_y_continuous(breaks=c(1,10,20,30))+
#                             #coord_cartesian(ylim=c(0,35),xlim=c(0.25,4.25))+
#                             scale_fill_manual(values = c("red" = "white", 'blue' = 'blue'), 
#                                        labels = c("red" = "More Glutamate Released", "blue" = "Saturated")) +
#                             scale_colour_manual(values = c("red" = "black", 'blue' = 'blue'), 
#                                          labels = c("red" = "More Glutamate Released", "blue" = "Saturated")) +                    
#                             theme_tufte()+
#                             my.theme+
#                             theme(legend.position = "none")#,
#                             #               legend.justification = c('right','top'),
#                             #               axis.text.x=element_text(colour="black", size=20, family="sans"),
#                             #               axis.text.y=element_text(colour="black", size=20, family="sans"),
#                             #               axis.title=element_text(colour="black", size=24, family="sans"),
#                             #               axis.ticks.length=unit(.25, "cm"),
#                             #               axis.ticks = element_line(size=1),
#                             #               strip.text = element_text(colour="black", size = 16, family = "sans")
#                             #               )

#             save_plot(filename = paste0(filename, "distribution_plot"), plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

#         }
            

# create_distribution_plots(vis_distributions, "Ca_mM", "mean_quanta", expression("["*Ca^'2+'*"] (mM)"), expression(E[Delta*"F/F"]/S[Delta*"F/F"]),"RRP_score", "Ca_vs_mean_quanta", "", FALSE,FALSE)

# create_distribution_plots(vis_distributions, "Ca_mM", "CV_amplitude", expression("["*Ca^'2+'*"] (mM)"), expression(CV[Delta*"F/F"]),"RRP_score", "Ca_vs_CV_amplitude", "", FALSE,FALSE)

# create_distribution_plots(vis_distributions, "Ca_mM", "releaseProbability_perROI", expression("["*Ca^'2+'*"] (mM)"), expression(P[iGlu]),"RRP_score", "Ca_vs_releaseProb", "", FALSE,FALSE)

# create_distribution_plots(vis_distributions, "Ca_mM", "var_amplitude", expression("["*Ca^'2+'*"] (mM)"), expression(sigma^2),"RRP_score", "Ca_vs_var_amplitude", "", FALSE,FALSE)

# create_distribution_plots(vis_distributions, "Ca_mM", "CV_tau_decay", expression("["*Ca^'2+'*"] (mM)"), expression(CV[tau*"_decay"]),"RRP_score", "Ca_vs_CV_tau_decay", "", FALSE,FALSE)

# create_distribution_plots(vis_distributions, "Ca_mM", "release_ratio", expression("["*Ca^'2+'*"] (mM)"), expression(frac(E/S["4Ca"], E/S["2Ca"])),"RRP_score", "Ca_vs_release_ratio", "", FALSE,FALSE)







# #     avgtrace_Plot<-ggplot(clean_peaks_normTime, aes(x=new_normTime, y = dFF,  group=Ca, colour=Ca))+
# #                                                     #geom_path(colour="black",size=0.4,alpha=0.2)+
# #                                                     stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=2.5,alpha=0.8)+
# #                                                     geom_vline(xintercept = 0, lty="longdash",colour='black',size=1)+
                                                    
# #                                                     labs( x="Normalized time (s)",
# #                                                           y=expression(Delta*"F/F"),
# #                                                           #title="Stimulus-evoked iGluSnFR3 activity",
# #                                                           colour=expression("[Ca"^{"2+"}*"], mM"),
# #                                                                         tag = "A")+
# #                                                     {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
# #                                                     {if(color_switch)scale_colour_manual(labels=new_levels,values=color_override)}+
# #                                                     scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
# #                                                     theme_tufte()+
# #                                                     my.theme+
# #                                                     #facet_grid(~Ca, labeller = label_parsed)+
# #                                                     theme(legend.title = element_text(colour="black", size=28, family="sans"),
# #                                                         legend.position=c(0.75,0.5))

# #                                                     #avg_peaks  <<- avg_PP_tracePlot
# #                                                      #       rm(avg_PP_tracePlot
                                                        

# #                                                      save_plot(filename=paste0("averagePeaks_basalRP"),plot=avgtrace_Plot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                         

# #                      #  averaged_traces <<- avgtrace_Plot
                         


#             #ggsave(filename=paste0("quanta_plot_v1"),plot=releaseProbPlot, device="jpeg",dpi=600, units="in",width=plot_height*1,height=plot_height*1)
            
# #            quanta_conv<<- releaseProbPlot


# ###generate for png plots

# corr_stats_saturation_score<- corr_stats %>% select(trackROI_key,RRP_score,release_ratio)


# RRP_score.labs <- c("Weaker Release","Stronger Release")
# names(RRP_score.labs) <- c("blue", "red")

# bar_data<- corr_stats %>% group_by(RRP_score) %>% summarise(n=n()) %>% mutate(freq=n/sum(n)) %>% mutate(asPercent = label_percent()(freq),
#                                                                                                         groupfact = "Boutons")


# #for printing
# check_bar_data<-bar_data %>% select(RRP_score,n,asPercent)
# print(check_bar_data)
# bar_data$RRP_score <- factor(bar_data$RRP_score,levels = c("red","blue"))
# #levels(bar_data$RRP_score) <- c('red','blue')

# ###REPRESENT SATURATERS VS ADDITIONAL RELEASERS AS FREQUENCY STACKED BAR CHART


# freq_bar_chart<-ggplot(bar_data, aes(x = groupfact, y=freq,fill = RRP_score)) +
#   geom_bar(position = "fill",stat="identity",colour="black",size=1,alpha=0.7) +
#   geom_label(aes(label = asPercent), 
#             stat = "identity",
#             position = position_fill(vjust = 0.5),
#             size=7,
#             fill="white",
#             show.legend=NA
#             )+
#   labs(x = expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "2 mM"~Ca^'2+')),
#         y = "Density",
#         tag = "") +
#   theme_tufte() +
#   coord_cartesian(ylim=c(0,1.1))+
#   scale_fill_manual(values = c("red" = "white", 'blue' = 'blue'), 
#                                         labels = c("red" = "More Glutamate Released", "blue" = "Saturated")) +
#   #guides()+
#   scale_y_continuous(expand=c(0,0),breaks=c(0,0.5,1.0))+
                      
#   my.theme +
#   theme(legend.position = "none",
#         axis.text.x=element_text(colour="white", size=scalefactor*34, family="sans"), #angle=45, hjust=1),
#         axis.title.x=element_text(colour="white", size=scalefactor*40, family="sans"),
#         axis.ticks.x=element_blank(),
#         axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank()
#         #axis.line=element_blank()                 
#         )

# save_plot(filename = "bar_test", plot=freq_bar_chart, plot_width=plot_height*2,plot_height=plot_height, scale_f = scalefactor,dpi_val = 300)
# bar_test<<- freq_bar_chart









# create_RRP_plot <- function(data, x_variable, y_variable, x_label,  y_label, group_var,filename,tag_var,add_text_labels,global_output) {
  

#             RRP_plot <- ggscatter(data, 
#                                x = x_variable, 
#                                y = y_variable,
#                                size = 3,
#                                color = group_var,
#                                shape = 21,
#                                stroke = 1.1,
#                                fill = group_var,
#                                alpha = 0.5) +
#                      {if(add_text_labels)geom_abline(lty = 'solid', intercept = 0, slope = 1)} +
#                      scale_fill_manual(values = c("red" = "white", 'blue' = 'blue'), 
#                                        labels = c("red" = "More Glutamate Released", "blue" = "Saturated")) +
#                      scale_colour_manual(values = c("red" = "black", 'blue' = 'blue'), 
#                                          labels = c("red" = "More Glutamate Released", "blue" = "Saturated")) +     
#                      {if(add_text_labels)geom_text(label = expression(atop("Response", "Saturation")), 
#                                y = 5, x = 22.5, size = 7, inherit.aes = FALSE, colour = 'blue')} +
#                      {if(add_text_labels)geom_text(label = expression(atop("Additional", "Glutamate Release")), 
#                                y = 30, x = 12.5, size = 7, inherit.aes = FALSE)} +
#                      coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
#                      labs(
#                        x = x_label,
#                        y = y_label,
#                        tag = tag_var
#                      ) +
#                      theme_tufte() +
#                      my.theme +
#                      theme(
#                        legend.position = "none"
#                      )
  
#   save_plot(
#     filename = paste0(filename, ""),
#     plot = RRP_plot,
#     plot_width = plot_height*1.5,
#     plot_height = plot_height,
#     scale_f = scalefactor,
#     dpi_val = 600
#   )
  
#   if(global_output == TRUE){
#     assign(paste0(filename, "_plot"), RRP_plot, envir = .GlobalEnv)}

#   return(RRP_plot)


#   }



#  #create_RRP_plot(corr_stats, "mean_quanta_1Ca", "mean_quanta_2Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "1 mM"~Ca^'2+')),expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "2 mM"~Ca^'2+')),"RRP_score", "RRP_2_vs_1Ca", "", FALSE)
#  #create_RRP_plot(corr_stats, "mean_quanta_1Ca", "mean_quanta_4Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "1 mM"~Ca^'2+')),expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),"RRP_score", "RRP_4_vs_1Ca", "", FALSE)
#  create_RRP_plot(corr_stats, "mean_quanta_2Ca", "mean_quanta_4Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "2 mM"~Ca^'2+')),expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),"RRP_score", "RRP_4Ca_vs_RRP_2Ca", "E", TRUE,TRUE)



# create_mean_variance_plot <- function(data, x_variable, y_variable, x_label,  y_label, group_var,filename,tag_var,add_text_labels,global_output) {



#                             max_y = 1.5


#                             # binomial model for uniform release probability
#                             # variance = IQ - I^2/ N 
#                             # y = Ax - Bx^2 where B is 1/N

#                             # binom_model <- function(mean_quanta, A, B, C){
#                             #     C + A*mean_quanta - ( (mean_quanta^2)/B )
#                             # }

#                             # if(str_detect(y_variable, "var_amplitude" ) == TRUE){
#                             #     max_y = 15
#                             #                                             model_data<- data %>% dplyr::filter(RRP_score == "blue")
#                             #                                             names(model_data)[names(model_data) == y_variable] <- 'var_amplitude'
#                             #                                             names(model_data)[names(model_data) == x_variable] <- 'mean_quanta'


#                             #                                             binomial_fit <- nls.multstart::nls_multstart(var_amplitude ~ binom_model(mean_quanta, A, B, C ),
#                             #                                                                                             data = model_data,
#                             #                                                                                             lower=c(A=0.3, B=0.1,C=0),                
#                             #                                                                                             upper=c(A=1.5, B=50,C=1),           
#                             #                                                                                             start_lower = c(A=0, B=0.1,C=0),
#                             #                                                                                             start_upper = c(A=1.5, B=50,C=1),
#                             #                                                                                             iter = 500,
#                             #                                                                                             supp_errors = "Y")

#                             #                                             print(summary(binomial_fit))
#                             #                                             predframe <- tibble(mean_quanta=seq(from=min(model_data$mean_quanta,na.rm=TRUE), to=max(model_data$mean_quanta,na.rm=TRUE), length.out = 1024)) %>%
#                             #                                                             mutate(var_amplitude = predict(binomial_fit, newdata = list(mean_quanta=.$mean_quanta)))
                                                                                      

#                             #                                                           check_nls_fit<- ggplot(model_data, aes(x=mean_quanta, y=var_amplitude)) +
#                             #                                                             geom_point(size=3) +
#                             #                                                             geom_line(data = predframe, aes(x=mean_quanta, y=var_amplitude))+
#                             #                                                             coord_cartesian(xlim=c(0,35),ylim=c(0,max_y))+
#                             #                                                             scale_y_continuous(breaks=c(0,0.5,1.0))+

#                             #                                                             #coord_flip()+
#                             #                                                             labs(x=x_label, #expression(E[Delta*"F/F"] / S[Delta*"F/F"]*" at 1 "*Ca^'2+'),
#                             #                                                                     y=y_label, #expression(CV[Delta*"F/F"]*" at 1 "*Ca^'2+'),
#                             #                                                                     tag=tag_var
#                             #                                                                     )+
#                             #                                                             theme_tufte()+
#                             #                                                             my.theme+
#                             #                                                             theme(legend.position = "none")
                                                              
                                                                                    

#                             #                                             #check_nls_object<- plot_nls(binomial_fit, model_data)
#                             #                                             save_plot(filename = paste0("check_nls_fit_",y_variable,"_xtra_release"),
#                             #                                                                                 plot = check_nls_fit,
#                             #                                                                                 plot_width = plot_height*1.5,
#                             #                                                                                 plot_height = plot_height,
#                             #                                                                                 scale_f = scalefactor,
#                             #                                                                                 dpi_val = 600)
#                             #                                         }
                            
#                             RRP_plot <- ggscatter(data, 
#                                             x = x_variable, y = y_variable, 
#                                             size = 3,
#                                             color=group_var,
#                                             shape=21,
#                                             stroke=1.1,
#                                             group=group_var,
#                                             fill=group_var,
#                                             alpha=0.5#,
#                                             #add = "reg.line", conf.int = TRUE,
#                                             #cor.coef = FALSE, cor.method = "pearson"
#                                            )+
#                             {if(add_text_labels)stat_regline_equation(label.x = 1, label.y=0.95*max_y,size=4,formula = y ~ poly(x, 2),
#                                                     aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))} +
#                             {if(add_text_labels)geom_smooth(method = "lm", formula = y ~ poly(x, 2))}+
#                             #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..r.label.., sep = "~")), label.x = 1, label.y=0.95*max_y, size=7, geom = "text")} +
#                             #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..p.label.., sep = "~")), label.x = 1, label.y=0.82*max_y, size=7, geom = "text")} +
                            
#                             #{if(add_text_labels)stat_regline_equation(label.x = 1, label.y=0.82*max_y,size=7)}+
                              
#                             scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
#                             scale_colour_manual(values = c("red"  = "black",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+     
#                             coord_cartesian(xlim=c(0,35),ylim=c(0,max_y))+
#                             scale_y_continuous(breaks=c(0,0.5,1.0))+

#                             #coord_flip()+
#                             labs(x=x_label, #expression(E[Delta*"F/F"] / S[Delta*"F/F"]*" at 1 "*Ca^'2+'),
#                                     y=y_label, #expression(CV[Delta*"F/F"]*" at 1 "*Ca^'2+'),
#                                     tag=tag_var
#                                     )+
#                             theme_tufte()+
#                             facet_grid(~RRP_score, labeller = labeller(RRP_score = RRP_score.labs))+
#                             my.theme+
#                             theme(legend.position = "none")
  
#   save_plot(
#     filename = paste0(filename, ""),
#     plot = RRP_plot,
#     plot_width = plot_height*1.5,
#     plot_height = plot_height,
#     scale_f = scalefactor,
#     dpi_val = 600
#   )

#   if(global_output == TRUE){
#     assign(paste0(filename, "_plot"), RRP_plot, envir = .GlobalEnv)}

#   return(RRP_plot)
# }

# create_mean_variance_plot(corr_stats, "mean_quanta_1Ca", "CV_amplitude_1Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "1 mM"~Ca^'2+')),expression(atop(CV[Delta*"F/F"], "1 mM"~Ca^'2+')),"RRP_score", "RRP_1_vs_CV_1Ca", "", TRUE,FALSE)
# create_mean_variance_plot(corr_stats, "mean_quanta_2Ca", "CV_amplitude_2Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "2 mM"~Ca^'2+')),expression(atop(CV[Delta*"F/F"], "2 mM"~Ca^'2+')),"RRP_score", "RRP_2_vs_CV_2Ca", "", TRUE,FALSE)
# create_mean_variance_plot(corr_stats, "mean_quanta_4Ca", "CV_amplitude_4Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),expression(atop(CV[Delta*"F/F"], "4 mM"~Ca^'2+')),"RRP_score", "RRP_4_vs_CV_4Ca", "", TRUE,FALSE)

# create_mean_variance_plot(corr_stats, "mean_quanta_4Ca", "CV_amplitude_1Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),expression(atop(CV[Delta*"F/F"], "1 mM"~Ca^'2+')),"RRP_score", "RRP_4_vs_CV_1Ca", "F", TRUE,TRUE)

# create_mean_variance_plot(corr_stats, "mean_amplitude_1Ca", "var_amplitude_1Ca", expression(atop(Delta*"F/F", "1 mM"~Ca^'2+')),expression(atop(sigma^2, "1 mM"~Ca^'2+')),"RRP_score", "RRP_1_vs_var_1Ca", "", TRUE,FALSE)
# create_mean_variance_plot(corr_stats, "mean_amplitude_2Ca", "var_amplitude_2Ca", expression(atop(Delta*"F/F", "2 mM"~Ca^'2+')),expression(atop(sigma^2, "2 mM"~Ca^'2+')),"RRP_score", "RRP_2_vs_var_2Ca", "", TRUE,FALSE)
# create_mean_variance_plot(corr_stats, "mean_amplitude_4Ca", "var_amplitude_4Ca", expression(atop(Delta*"F/F", "4 mM"~Ca^'2+')),expression(atop(sigma^2, "4 mM"~Ca^'2+')),"RRP_score", "RRP_4_vs_var_4Ca", "", TRUE,FALSE)



# create_violin_pairs_plot <- function(data, y_variable, y_label, group_var,filename,tag_var,add_text_labels,global_output) {



#                     max_y = 35

#                      plot <- ggplot(data, aes_string(x = group_var, y = y_variable, group = group_var, colour = group_var)) +
#                                 geom_sina(size=1.5,alpha=0.5)+
#                                 geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=1) +
                                                

#                                 labs(
#                                   y = y_label,
#                                   x = "",
#                                   colour = "",
#                                   tag = tag_var
#                                 ) +
#                                 #scale_y_continuous(breaks = c(0, 0.5, 1.0,1.5)) + #expand=c(0,0)
#                                 #coord_cartesian(ylim = c(0, #y_upper_bound)) +
#                                 #scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
#                                 scale_colour_manual(values = c("red"  = "black",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+     
#                                 theme_tufte() +
#                                 my.theme +
#                                 theme(legend.position = "none")



#   save_plot(
#     filename = paste0(filename, ""),
#     plot = plot,
#     plot_width = plot_height*1.5,
#     plot_height = plot_height,
#     scale_f = scalefactor,
#     dpi_val = 600
#   )

#   if(global_output == TRUE){
#     assign(paste0(filename, "_plot"), plot, envir = .GlobalEnv)}

#   return(plot) 
       
# }

# create_violin_pairs_plot(corr_stats, "CV_amplitude_1Ca", expression(atop(CV[Delta*"F/F"], "1 mM"~Ca^'2+')),"RRP_score", "RRP_violins_vs_CV_1Ca", "", TRUE,FALSE)
# create_violin_pairs_plot(corr_stats, "CV_amplitude_2Ca", expression(atop(CV[Delta*"F/F"], "2 mM"~Ca^'2+')),"RRP_score", "RRP_violins_vs_CV_2Ca", "", TRUE,FALSE)
# create_violin_pairs_plot(corr_stats, "CV_amplitude_4Ca", expression(atop(CV[Delta*"F/F"], "4 mM"~Ca^'2+')),"RRP_score", "RRP_violins_vs_CV_4Ca", "", TRUE,FALSE)
# create_violin_pairs_plot(corr_stats, "releaseProbability_perROI", expression(atop(P[iGlu], "0.5 mM"~Ca^'2+')),"RRP_score", "RRP_violins_vs_PiGlu_05Ca", "", TRUE,FALSE)
# create_violin_pairs_plot(corr_stats, "mean_quanta_4Ca", expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),"RRP_score", "RRP_violins_vs_mean_quanta_4Ca", "", TRUE,FALSE)









# max_y = 1.1

#                     releaseProbPlot<-ggscatter(corr_stats, y= "releaseProbability_perROI", x = "mean_quanta_4Ca",  
#                                             size = 3,
#                                             color='RRP_score',
#                                             shape=21,
#                                             stroke=1.1,
#                                             group="RRP_score",
#                                             fill="RRP_score",
#                                             alpha=0.5, 
#                                             add = "reg.line", conf.int = TRUE,
#                                             #add.params = list(color='blue',fill='lightgrey'), 
#                                             cor.coef = FALSE, cor.method = "pearson"#,
#                                             )+
#                             ggpubr::stat_cor(aes(label = paste(..r.label.., sep = "~")), label.x = 0.1, label.y=0.95*max_y, size=7, geom = "text") +
#                             ggpubr::stat_cor(aes(label = paste(..p.label.., sep = "~")), label.x = 0.1, label.y=0.82*max_y, size=7, geom = "text") +
#                             scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
#                             scale_colour_manual(values = c("red"  = "black",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+     
                            
#                             labs(   x=expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),
#                                      y=expression(atop(P[iGlu], "0.5 mM"~Ca^'2+')),
                                  
#                                    # title="",
#                                     tag="G"
#                                     )+
#                             coord_cartesian(xlim=c(0,35),ylim=c(0,max_y))+
#                             scale_y_continuous(breaks=c(0,0.5,1))+
#                             theme_tufte()+
#                             facet_grid(~RRP_score, labeller = labeller(RRP_score = RRP_score.labs))+
#                             my.theme+
#                             theme(legend.position = "none")
                           

#             save_plot(filename = "RRP_vs_lowCa_prob", plot=releaseProbPlot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)



#                 releaseProbPlot<-ggscatter(corr_stats, y= "CV_amplitude_1Ca", x = "mean_quanta_4Ca",  
#                                             size = 3,
#                                             color='RRP_score',
#                                             shape=21,
#                                             stroke=1.1,
#                                             group="RRP_score",
#                                             fill="RRP_score",
#                                             alpha=0.5, 
#                                             add = "reg.line", conf.int = TRUE,
#                                             #add.params = list(color='blue',fill='lightgrey'), 
#                                             cor.coef = FALSE, cor.method = "pearson"#,
#                                             )+
#                             ggpubr::stat_cor(aes(label = paste(..r.label.., sep = "~")), label.x = 0.1, label.y=0.95*max_y, size=7, geom = "text") +
#                             ggpubr::stat_cor(aes(label = paste(..p.label.., sep = "~")), label.x = 0.1, label.y=0.82*max_y, size=7, geom = "text") +
#                             scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
#                             scale_colour_manual(values = c("red"  = "black",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+     
                            
#                             labs(   x=expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4 mM"~Ca^'2+')),
#                                      y=expression(atop(CV[Delta*"F/F"], "1 mM"~Ca^'2+')),
                                  
#                                    # title="",
#                                     tag="G"
#                                     )+
#                             coord_cartesian(xlim=c(0,35),ylim=c(0,max_y))+
#                             scale_y_continuous(breaks=c(0,0.5,1))+
#                             theme_tufte()+
#                             facet_grid(~RRP_score, labeller = labeller(RRP_score = RRP_score.labs))+
#                             my.theme+
#                             theme(legend.position = "none")
                           

#             save_plot(filename = "RRP_vs_CV_1Ca_prob", plot=releaseProbPlot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
#             RRP_4Ca_vs_PiGlu_05<<- releaseProbPlot

#             rm(max_y)



#####generate releaseProb_calc version of N_sites

# library(nls.multstart) 

# #check<- releaseProb_calc %>% ungroup() %>% select(trackROI_key,mean_amplitude, var_amplitude) %>% dplyr::filter(trackROI_key == "GluSnFR3-dish05-plate09-region02-ROI0006")
# #print(check)
# active_4Ca_list<- releaseProb_calc %>% dplyr::filter(Ca_mM == 4, releaseProbability_perROI == 1, !is.na(var_amplitude) )
# active_4Ca_list<- unique(active_4Ca_list$trackROI_key)

# active_2Ca_list<- releaseProb_calc %>% dplyr::filter(Ca_mM == 2, releaseProbability_perROI >0.1, !is.na(var_amplitude) )
# active_2Ca_list<- unique(active_2Ca_list$trackROI_key)

# active_1Ca_list<- releaseProb_calc %>% dplyr::filter(Ca_mM == 1, releaseProbability_perROI > 0.1, !is.na(var_amplitude) )
# active_1Ca_list<- unique(active_1Ca_list$trackROI_key)


# print(length(active_1Ca_list))
# #### nonuniform parameters 
# #
# # y = Q* x - (Q*x^2(1+alpha)/(x + N*Q*alpha))

# # here we can solve directly for N... ?
# # Leave q as a free parameter to see what it comes up with
# # alpha is a free parameter

# # linear_reg_model<- releaseProb_calc %>% ungroup() %>% group_by(trackROI_key) %>% dplyr::filter(Ca_mM != 0.5) %>%
# #                         select(trackROI_key, mean_amplitude,var_amplitude) %>% 
# #                         dplyr::filter(!is.na(var_amplitude), !is.na(mean_amplitude), trackROI_key %in% active_2Ca_list) %>%
# #                         nest() %>%
# #                          mutate(binom_fit = purrr::map(data, ~nls_multstart(var_amplitude ~ ( ( Q*mean_amplitude - (Q*mean_amplitude^2 * (1 + alpha) )/(mean_amplitude + N*Q*alpha) )*(1+0.79^2)) ,
# #                                                         data=.x,
# #                                                         iter = 500,
# #                                                         start_lower = c(Q = 0.1, N = 1, alpha = 0),
# #                                                         start_upper = c(Q = 1.5, N = 100, alpha = 1000),
# #                                                         supp_errors = 'Y',
# #                                                         na.action = na.omit,
# #                                                         lower = c(Q = 0, N = 0, alpha = 0))),
# #                                 tidied = map(binom_fit, tidy),
# #                                 augment = map(binom_fit, augment)
# #                                 ) #%>%

# N_sites_nonunif_reg_model<- releaseProb_calc %>% ungroup() %>% group_by(trackROI_key) %>% dplyr::filter(Ca_mM != 0.5) %>%
#                         select(trackROI_key, mean_amplitude,var_amplitude) %>% 
#                         #dplyr::filter(!is.na(var_amplitude), !is.na(mean_amplitude), trackROI_key %in% active_4Ca_list, trackROI_key %in% active_2Ca_list, trackROI_key %in% active_1Ca_list) %>%
#                         nest() %>%
#                          mutate(binom_fit = purrr::map(data, ~nls_multstart(var_amplitude ~ ( ( Q*mean_amplitude - (Q*mean_amplitude^2 * (1 + alpha) )/(mean_amplitude + N*Q*alpha) )*(1+0.79^2)) ,
#                                                         data=.x,
#                                                         iter = 500,
#                                                         start_lower = c(Q = 0.1, N = 1, alpha = 0),
#                                                         start_upper = c(Q = 1.5, N = 100, alpha = 1000),
#                                                         supp_errors = 'Y',
#                                                         na.action = na.omit,
#                                                         lower = c(Q = 0, N = 0, alpha = 0))),
#                                 tidied = map(binom_fit, tidy),
#                                 augment = map(binom_fit, augment)
#                                 ) #%>%
#                         #unnest(tidied) %>%
#                         #select(trackROI_key,term,estimate) %>%
#                         #spread(term,estimate) %>%

# check_nonunif_model <- N_sites_nonunif_reg_model %>% unnest(tidied) %>%
#                         select(trackROI_key,term,estimate) %>%
#                         spread(term,estimate) #%>%
#                         #mutate(Pr = mean_amplitude/(N*Q))
                        

# generate_nonunif_fit_predict<- check_nonunif_model %>% group_by(across(everything())) %>%
#                                               summarise(amp_points = seq(from = 0, to = 15, length.out=150)) %>%
#                                               ungroup() %>%
#                                               group_by(trackROI_key) %>% 
#                                               mutate(var_pred = ( Q*amp_points - (Q*amp_points^2 * (1 + alpha) )/(amp_points + N*Q*alpha) )*(1+0.79^2) ) 

# names(generate_nonunif_fit_predict)[names(generate_nonunif_fit_predict) == 'amp_points'] <- 'mean_amplitude'
# names(generate_nonunif_fit_predict)[names(generate_nonunif_fit_predict) == 'var_pred'] <- 'var_amplitude'



########################################################################

#trying to produce an offset to correct mean_variance plots.
#mean_amplitude_06Ca subtracted from all other ampltiudes for centering mean/variance at 0

##########################################################################

offset_06Ca<- corr_stats %>% select(trackROI_key, mean_amplitude_06Ca)

corrected_releaseProb_calc<- left_join(releaseProb_calc,offset_06Ca) %>% mutate(normalized_mean_amplitude = mean_amplitude - mean_amplitude_06Ca)









                                                            ## jail
adj_nan_releaseProb_calc<- corrected_releaseProb_calc %>% dplyr::filter(!trackROI_key %in% c("GluSnFR3-dish02-plate08-region01-ROI0016", "GluSnFR3-dish03-plate02-region01-ROI0015","GluSnFR3-dish03-plate08-region01-ROI0006"))
adj_nan_releaseProb_calc[is.na(adj_nan_releaseProb_calc)] <- 0
adj_nan_releaseProb_calc[is.nan(adj_nan_releaseProb_calc)] <- 0

N_sites_reg_model<- adj_nan_releaseProb_calc %>% ungroup() %>% group_by(trackROI_key) %>% dplyr::filter(Ca_mM != 0.5) %>%
                        select(trackROI_key, mean_amplitude,var_amplitude) %>% 
                       # dplyr::filter(trackROI_key %in% active_4Ca_list) %>%
                        nest() %>%
                         mutate(binom_fit = purrr::map(data, ~nls_multstart(var_amplitude ~ corr_factor + A*mean_amplitude - B*mean_amplitude^2,
                                                        data=.x,
                                                        iter = 1000,
                                                        start_lower = c(A = 0.001, B = 0.001, corr_factor = -20),
                                                        start_upper = c(A = 1000, B = 1000, corr_factor = 20),
                                                        supp_errors = 'Y',
                                                        na.action = na.omit,
                                                        lower = c(A = 0, B= 0, corr_factor = -20))),
                                tidied = map(binom_fit, tidy),
                                augment = map(binom_fit, augment)
                                ) #%>%
                        #unnest(tidied) %>%
                        #select(trackROI_key,term,estimate) %>%
                        #spread(term,estimate) %>%
                        #mutate(N_min = 1/B)
                    

check_good_N_sites <- N_sites_reg_model %>% unnest(tidied) %>%
                        select(trackROI_key,term,estimate) %>%
                        spread(term,estimate) %>%
                        mutate(N_min = 1/B) %>%
                    #dplyr::filter( A > 0 ,B > 0, N_min < 10) %>%
                    mutate(N_sites_group = ifelse(N_min < 20, "A", "B"))

print(check_good_N_sites[,c(1:5)])

good_N_sites<- adj_nan_releaseProb_calc %>% dplyr::filter(trackROI_key %in% unique(check_good_N_sites$trackROI_key))
good_N_sites <- left_join(good_N_sites, check_good_N_sites) 


generate_fit_predict<- check_good_N_sites %>% group_by(across(everything())) %>%
                                              summarise(amp_points = seq(from = 0, to = 15, length.out=150)) %>%
                                              ungroup() %>%
                                              group_by(trackROI_key) %>% 
                                              mutate(var_pred = corr_factor + A*amp_points - B*amp_points^2,
                                                     N_sites_group = ifelse(N_min < 20, "A", "B"))




generate_fit_predict <- generate_fit_predict %>% group_by(trackROI_key) %>% select(trackROI_key,  amp_points, var_pred) #N_sites_group,
names(generate_fit_predict)[names(generate_fit_predict) == 'amp_points'] <- 'mean_amplitude'
names(generate_fit_predict)[names(generate_fit_predict) == 'var_pred'] <- 'var_amplitude'

check_fitted<-N_sites_reg_model  %>% unnest(augment)
get_sum_sq_resid <- check_fitted %>% group_by(trackROI_key) %>%
                                  summarise( sum_sq_resid = sum(.resid^2,na.rm=TRUE))

keep_ROIs  <- get_sum_sq_resid #%>% dplyr::filter(!is.na(sum_sq_resid) ) 



print("info on residuals: ")
print("")
print("Max")
print(max(get_sum_sq_resid$sum_sq_resid, na.rm=TRUE))
print("Min")
print(min(get_sum_sq_resid$sum_sq_resid, na.rm=TRUE))

print("Median")
print(median(get_sum_sq_resid$sum_sq_resid, na.rm=TRUE)) 

print("Number of ROIs retained after filtering on sum_sq_resid < 1 ")
print(length(unique(keep_ROIs$trackROI_key)))

print(check_fitted)


print(N_sites_reg_model)

good_N_sites <- good_N_sites %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key))
generate_fit_predict <- generate_fit_predict %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key))



y_limit = 1.5*max(unique(releaseProb_calc$var_amplitude),na.rm=TRUE)






N_sites<-ggplot(generate_fit_predict, aes(x=mean_amplitude, y=var_amplitude, group=trackROI_key))+ 
                                            #x = "mean_amplitude", y = "var_amplitude", 
                                            #size = 2,
                                            #color='Ca',
                                            #group='trackROI_key',
                                            #shape=21,
                                            #stroke=1.1,
                                            #fill='white',
                                            #alpha=0.5
                                            #)+
                            geom_line( colour="blue", size=1, alpha=0.5)+
                            geom_point(data=good_N_sites,aes(x=mean_amplitude,y=var_amplitude, colour=Ca), size=2, shape=21, stroke=1.1,fill=NA,alpha=0.5)+ 
                            
                            #geom_vline(xintercept = 1, lty='dashed',size=0.5,colour="black")+
                            #geom_smooth(aes(group=trackROI_key),method = "lm", formula = y ~ poly(x, 2), se = FALSE,alpha=0.1)+
                            
                           # stat_regline_equation(label.x = 1, label.y=0.95*y_limit,size=7,formula = y ~ poly(x, 2),
                            #                        aes(label =  paste(..eq.label.., sep = "~~~~"))) +
                            #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..r.label.., sep = "~")), label.x = 1, label.y=0.95*max_y, size=7, geom = "text")} +
                            #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..p.label.., sep = "~")), label.x = 1, label.y=0.82*max_y, size=7, geom = "text")} +
                            
                            #{if(add_text_labels)stat_regline_equation(label.x = 1, label.y=0.82*max_y,size=7)}+
                              
                            #scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
                            scale_colour_manual(values = color_override)+
                            scale_x_continuous(expand=c(0,0))+
                            scale_y_continuous(expand=c(0,0))+     
                            coord_cartesian(xlim=c(0,15),ylim=c(0,y_limit))+
                            #scale_y_continuous(breaks=c(0,0.5,1.0))+

                            #coord_flip()+
                            labs(x=expression("Mean"~Delta*"F/F"),
                                    y=expression("Variance"~Delta* "F/F"^2),
                                    tag="",
                                    title="Population Mean-Variance"
                                    )+
                            theme_tufte()+
                            #facet_grid(~N_sites_group)+
                            my.theme+
                            theme(legend.position = "none")

                            save_plot(filename = "mean-variance_nlsFits", plot=N_sites, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 300)
            





subDir <- "check_N_sites"   


source(paste0(path,"sub-functions/setFigureDirectory.R") ) 

good_ROIs<- unique(keep_ROIs$trackROI_key)

for(i in 1:length(good_ROIs)){



        current_ROI = good_ROIs[i]

        tmp_fit_predict<- generate_fit_predict %>% dplyr::filter(trackROI_key == current_ROI)
        tmp_N_sites <- good_N_sites %>% dplyr::filter(trackROI_key == current_ROI)
        tmp_fit_values<- check_good_N_sites %>% dplyr::filter(trackROI_key == current_ROI)

        # tmp_nonunif_fit_predict<- generate_nonunif_fit_predict %>% dplyr::filter(trackROI_key == current_ROI)
        # tmp_nonunif_fit_values<- check_nonunif_model %>% dplyr::filter(trackROI_key == current_ROI)


       # if(unique(tmp_nonunif_fit_values$N) > 1000){
         #   add_second_model = FALSE
        #} else {
            add_second_model = FALSE
        #}

        x_limit = 2*max(unique(tmp_fit_predict$mean_amplitude),na.rm=TRUE)
        y_limit = 1.5*max(unique(tmp_fit_predict$var_amplitude),na.rm=TRUE)


                    N_sites<-ggplot(tmp_fit_predict, aes(x=mean_amplitude, y=var_amplitude, group=trackROI_key))+ 
                                                                
                                                geom_line( colour="black", size=1, alpha=0.6, lty="longdash")+
                                                {if(add_second_model)geom_line(data=tmp_nonunif_fit_predict, aes(x=mean_amplitude,y=var_amplitude, group=trackROI_key), colour="blue", size=1, lty='dotdash',alpha=0.7)}+
                                                #geom_smooth(data=tmp_N_sites, aes(x=mean_amplitude, y=var_amplitude, group=trackROI_key),method = "lm", formula = y ~ x, se = FALSE,alpha=0.5,colour="forestgreen")+
                                                # stat_regline_equation(label.x = x_limit*0.6, label.y=0.67*y_limit,size=7,formula = y ~ x, colour='forestgreen',
                                                 #                       aes(label =  paste(..eq.label.., sep = "~~~~"))) +
                                               
                                                geom_smooth(data=tmp_N_sites, aes(x=mean_amplitude, y=var_amplitude, group=trackROI_key),method = "lm", formula = y ~ poly(x, 2), se = FALSE,alpha=0.5,colour="red")+
                                                stat_regline_equation(label.x = x_limit*0.6, label.y=0.60*y_limit,size=7,formula = y ~ poly(x, 2), colour='red',
                                                                        aes(label =  paste(..eq.label.., sep = "~~~~"))) +
                                
                                                geom_point(data=tmp_N_sites,aes(x=mean_amplitude,y=var_amplitude, colour=Ca), size=3, shape=21, stroke=1.5,fill="white",alpha=1)+
                                                
                                                
                                                #unif
                                                geom_text(label = paste0("A = ", round(unique(tmp_fit_values$A),2 )), size=7, x= x_limit*0.7, y=y_limit*0.95) +
                                                geom_text(label = paste0("B = ", round(unique(tmp_fit_values$B),2 )), size=7, x= x_limit*0.7, y=y_limit*0.87) +
                                                geom_text(label = paste0("N_sites = ", round(unique(tmp_fit_values$N_min),2 )), size=7, x= x_limit*0.7, y=y_limit*0.80) +
                                                
                                                #nonunif
                                                {if(add_second_model)geom_text(label = paste0("Q = ", round(unique(tmp_nonunif_fit_values$Q),2 )), size=7, x= x_limit*0.7, y=y_limit*0.53, colour='blue')} +
                                                #{if(add_second_model)geom_text(label = paste0("alpha = ", round(unique(tmp_nonunif_fit_values$alpha),2 )), size=7, x= x_limit*0.7, y=y_limit*0.53,colour='blue')} +
                                                {if(add_second_model)geom_text(label = paste0("N_sites = ", round(unique(tmp_nonunif_fit_values$N),2 )), size=7, x= x_limit*0.7, y=y_limit*0.46,colour='blue')} +
                                                

                                                
                                             
                                                scale_colour_manual(values = color_override)+
                                                scale_x_continuous(expand=c(0,0),breaks=c(0,2,5,10,15,20,30))+
                                                scale_y_continuous(expand=c(0,0))+     
                                                coord_cartesian(xlim=c(0,x_limit),ylim=c(0,y_limit))+
                                                #scale_y_continuous(breaks=c(0,0.5,1.0))+

                                                #coord_flip()+
                                                labs(x=expression(Delta*"F/F"),
                                                        y=expression(sigma[Delta*"F/F"]^2),
                                                        tag="",
                                                        title="Mean-Variance plot with binomial fit"#,
                                                        #subtitle=current_ROI
                                                        )+
                                                theme_tufte()+
                                                my.theme+
                                                facet_grid(~trackROI_key)+
                                                theme(legend.position = "none",
                                                        plot.margin = unit(c(1,1,1,1),"cm"))

                                                save_plot(filename = paste0("mean-variance_nlsFits_", current_ROI), plot=N_sites, plot_width=plot_height*2,plot_height=plot_height*1.2, scale_f = scalefactor,dpi_val = 300)
}


                            #nam <- paste("tmp_N_sites", count, sep = "")
                            #assign(nam, tmp_N_sites,envir = .GlobalEnv)
            
# ###### 1) PLOT RELEASE PROBABILITY PER ROI  #######

                    wrap_val = 50
                    yval='releaseProbability_perROI'

                    x_lower_bound = 0 #min(BasalRP_calc$Ca_mM)
                    x_upper_bound = max(releaseProb_calc$Ca_mM) *1.5
                    x_range = x_upper_bound - x_lower_bound
                    dot_y_range = max(releaseProb_calc$releaseProbability_perROI)/20
                    jitter_x_range = x_upper_bound/20

                    omit_inactive<- releaseProb_calc %>% dplyr::filter(Ca == "8Ca", releaseProbability_perROI == 0 )
                    omit_inactive_list<- unique(omit_inactive$trackROI_key)

                    fix_releaseProb_calc<- releaseProb_calc %>% dplyr::filter(!trackROI_key %in% omit_inactive_list)

                    fix_releaseProb_calc$groupfact <- factor(fix_releaseProb_calc$Ca_mM)    
                    fix_releaseProb_calc$Ca<- as.factor(as.character( (fix_releaseProb_calc$Ca) ) )
                    levels(fix_releaseProb_calc$Ca) <- new_levels

                                            

                    releaseProbPlot<-ggplot(fix_releaseProb_calc, aes(x=Ca_mM, y=releaseProbability_perROI))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=1.2,alpha=0.35)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            
                            #geom_smooth(formula=y~x, method='loess',alpha=1,se=FALSE,size=3,colour="black")+              
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white", alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression(P[iGlu]),
                                   # title="",
                                    tag="A"
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                            #guides(colour="none")+
                            scale_x_continuous(breaks=c(0.6,1,2,4,8),labels=c("0.6","1","2","4","8"))+
                            scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
                            coord_cartesian(ylim=c(0,1.0),xlim=c(0.25,8.25))+
                            guides(colour = guide_legend(override.aes = list(colour = color_override,
                                                                                size = c(4,4,4,4),
                                                                                alpha=c(1,1,1,1)
                                                                                )
                                                                            )
                                                                        )+                    
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "right")#,
                                    #legend.text=element_text(colour="black", size=scalefactor*18, family="sans"))#,
                            #               legend.justification = c('right','top'),
                            #               axis.text.x=element_text(colour="black", size=20, family="sans"),
                            #               axis.text.y=element_text(colour="black", size=20, family="sans"),
                            #               axis.title=element_text(colour="black", size=24, family="sans"),
                            #               axis.ticks.length=unit(.25, "cm"),
                            #               axis.ticks = element_line(size=1),
                            #               strip.text = element_text(colour="black", size = 16, family = "sans")
                            #               )


            #ggsave(filename=paste0("releaseProbPlot_v1"),plot=releaseProbPlot, device="jpeg",dpi=600, units="in",width=plot_height*1,height=plot_height*1)
            save_plot(filename = "releaseProbPlot_v1", plot=releaseProbPlot, plot_width=plot_height*2,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                peak_probability<<- releaseProbPlot



                releaseProbPlot<-ggplot(fix_releaseProb_calc, aes(x=Ca_mM, y=mean_quanta))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=1.2,alpha=0.35)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            
                            geom_hline(yintercept = 1, size=0.7,lty="dashed",colour="black")+
                            #geom_smooth(formula=y~x, method='loess',alpha=1,se=FALSE,size=3,colour="black")+              
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white",,alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression(E[Delta*"F/F"] / S[Delta*"F/F"]),
                                    tag="B"
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                            #guides(colour="none")+
                            scale_x_continuous(breaks=c(0.6,1,2,4,8),labels=c("0.6","1","2","4","8"))+
                            scale_y_continuous(breaks=c(1,5,10,15,20,30))+
                            coord_cartesian(ylim=c(0,17),xlim=c(0.25,8.25))+
                                                
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")#,
                            #               legend.justification = c('right','top'),
                            #               axis.text.x=element_text(colour="black", size=20, family="sans"),
                            #               axis.text.y=element_text(colour="black", size=20, family="sans"),
                            #               axis.title=element_text(colour="black", size=24, family="sans"),
                            #               axis.ticks.length=unit(.25, "cm"),
                            #               axis.ticks = element_line(size=1),
                            #               strip.text = element_text(colour="black", size = 16, family = "sans")
                            #               )

            save_plot(filename = "quanta_plot_v1", plot=releaseProbPlot, plot_width=plot_height*2,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
            
            #ggsave(filename=paste0("quanta_plot_v1"),plot=releaseProbPlot, device="jpeg",dpi=600, units="in",width=plot_height*1,height=plot_height*1)
            
            quanta_conv<<- releaseProbPlot
             





# #"Very~active:~P[iGlu]~>~0.75"
# #"Active:~0.75~>~P[iGlu]~>~0.25",
# #"Rarely~active:~0.25~>~P[iGlu]~>~0",
# #"Inactive:~P[iGlu]~=~0"),

# # ###### 7) PLOT ACTIVE/INACTIVE SYNAPSE PROPORTIONS  #######
# cc <- scales::seq_gradient_pal("blue", "orange", "Lab")(seq(0,1,length.out=6))
# parse.labels <- function(x) parse(text = x)
# blank_theme <- theme_minimal()+
#   theme(
#   axis.title.x = element_blank(),
#   axis.title.y = element_blank(),
#   panel.border = element_blank(),
#   panel.grid=element_blank(),
#   axis.ticks = element_blank(),
#   plot.title=element_blank()#element_text(size=14, face="bold")
#   )


# pie_chart_data<- releaseProb_calc %>% ungroup() %>% group_by(Ca, Ca_mM, groupfact,activity) %>% summarise(counts = length(trackROI_key))
# pie_chart_total<- pie_chart_data %>% group_by(Ca,Ca_mM,groupfact) %>% summarise(total_counts = sum(counts))

# pie_chart_final<- left_join(pie_chart_data, pie_chart_total) %>% mutate(proportion = counts/total_counts,
#                                                                         asPercent = proportion*100)

# print(pie_chart_final)

# pie_chart_final$activity = factor(pie_chart_data$activity, levels = c("inactive","rarely active","somewhat active","active","very active","multivesicular"))#("Inactive:~P[iGlu]~=~0","Rarely~active:~0.25~>~P[iGlu]~>~0", "Active:~0.75~>~P[iGlu]~>~0.25","Very~active:~P[iGlu]~>~0.75"))

# activity_labels = c(expression("Inactive: "*P[iGlu]~"="~0),expression("Rarely Active: "*0~"<"~P[iGlu]~"<"~0.25),expression("Somewhat Active: "*0.25~"\u2264"~P[iGlu]~"<"~0.50),expression("Active: "*0.50~"\u2264"~P[iGlu]~"<"~0.75),expression("Very Active: "*0.75~"\u2264"~P[iGlu]),expression("Multivesicular") )

# levels(pie_chart_final$activity) = activity_labels

# pie<- ggplot(pie_chart_final, aes(x="", y=asPercent, fill=activity))+
#         geom_bar(width = 1, stat = "identity",colour="black",size=0.05)+
#         coord_polar("y", start=0)+
#         scale_fill_manual(values=cc,labels=activity_labels)+
#         #labs(tag="B")+
#         guides(fill = guide_legend(reverse = TRUE,label.position="left"))+
#         facet_wrap(~Ca, ncol=2, labeller = label_parsed )+
#         #my.theme+
#         blank_theme +
#         theme(plot.tag = element_text(colour="black", size=scalefactor*48, family="sans",face="bold"),
#                 legend.position="left",
#                 legend.justification = c("left","center"),
#                 axis.text = element_blank(),
#                 legend.title=element_blank())


#             # pirid <- ggplot(pie, aes(x = "", y = Counts, fill = Puncta)) +
#             #           geom_col(color = "black") +
#             #           geom_text(aes(label = percent),
#             #                     position = position_stack(vjust = 0.5)) +
#             #           coord_polar(theta = "y") + theme(axis.ticks = element_blank(),
#             #                 axis.title = element_blank(),
#             #                 axis.text = element_text(size = 15), 
#             #                 legend.position = "none", # Removes the legend
#             #                 panel.background = element_rect(fill = "white")) + 
#             #           scale_color_viridis(option = "D", discrete = TRUE) +
#             #           #scale_fill_viridis(discrete = TRUE) +
#             #           theme_void()#+
#             #           #my.theme()
# save_plot(filename=paste0("piechart_trial"),plot=pie, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

# #pie_activity<<- pie
            




#                         check_inactives <- releaseProb_calc %>% dplyr::filter(activity == "inactive", Ca_mM == 4) 
#                         print("Inactive synapses were found in these imaging regions at Ca = 4 mM.")
#                         print(unique(check_inactives$imaging_region) )
#                         print("Inactive synapses were found at these ROIs at Ca = 4 mM.")
#                         print(unique(check_inactives$trackROI_key) )

#                         check_actives <- releaseProb_calc %>% dplyr::filter(activity == "very active", Ca_mM == 0.5) 
#                         print("Very active synapses were found in these imaging regions at Ca = 0.5 mM.")
#                         print(unique(check_inactives$imaging_region) )
#                         print("Very active synapses were found at these ROIs at Ca = 0.5 mM.")
#                         print(unique(check_inactives$trackROI_key) )
                            

#                         BasalRP_summary <- suppressMessages(releaseProb_calc %>% group_by(Ca,Ca_mM) %>% summarise(total_neuronCount = n_distinct(imaging_region),
#                                                                                                                     totalSynapseCount = n_distinct(trackROI_key, na.rm=TRUE),
#                                                                                                                     releaseProb = mean(releaseProbability_perROI),
#                                                                                                                     releaseProb_active = mean(releaseProbability_perROI[which(releaseProbability_perROI!=0)]),
#                                                                                                                     inactiveSynapseCount = length(releaseProbability_perROI[which(releaseProbability_perROI<=0.05)]),
#                                                                                                                     inactiveSynapseProportion = inactiveSynapseCount/totalSynapseCount,
#                                                                                                                     activeSynapseCount = length(releaseProbability_perROI[which(releaseProbability_perROI>0.05)]),
#                                                                                                                     activeSynapseProportion = activeSynapseCount/totalSynapseCount#,
#                                                                                                                    # mean_time_to_peak = mean(mean_time_to_peak, na.rm=TRUE)
#                                                                                                                )
#                                                             )   
                         
#                         print(paste0("This is the total number of neurons imaged in the dataset probably: ", BasalRP_summary$total_neuronCount) )
#                         print(paste0("This is the total number of synapses in the dataset probably: ", BasalRP_summary$totalSynapseCount) )
#                         print(paste0("Total inactive synapses is :", BasalRP_summary$inactiveSynapseCount))


#                          synapseProp_pivot<- tidyr::pivot_longer(BasalRP_summary, cols = c('inactiveSynapseProportion','activeSynapseProportion'), names_to = 'synapseProp', values_to ="value")

#                          synapseProp_pivot$synapseProp<- factor(synapseProp_pivot$synapseProp, levels = c('activeSynapseProportion','inactiveSynapseProportion'))
#                          synapse_labels = c("Active Boutons","Inactive Boutons")
                        
#                              plotval = 'Ca_mM'
#                              yval = 'value'
#                              fillval = 'synapseProp'

#                              # barPlot <-ggplot(synapseProp_pivot, aes_string(x=plotval,y=yval,colour=fillval)) +
#                              #                    geom_hline(yintercept = 0.1, size=0.7,lty="dashed",colour="black")+
                            
#                              #                    geom_line(alpha=1,size=1.5,lty="solid")+              
                            
#                              #                    #geom_line(size=2)+ 
#                              #                    geom_point(shape=21,stroke=1.1,size=4,fill="white")+
                                                
#                              #                    labs( x=expression("["*Ca^'2+'*"] (mM)"),
#                              #                          y=str_wrap("Proportion of boutons", wrap_val),
#                              #                          tag = "B")+
#                              #                    scale_colour_manual(labels = synapse_labels, values=c('forestgreen','grey24'))+
#                              #                    scale_x_continuous(breaks=c(0.5,1,2,4))+
#                              #                    scale_y_continuous(breaks=c(0,0.1,0.25,0.5,0.75,1.0),labels = number_format(accuracy = 0.01))+
#                              #                    coord_cartesian(ylim=c(0,1))+
#                              #                    theme_tufte()+
#                              #                    my.theme+
#                              #                    theme(legend.spacing.y = unit(2.0,'cm'),
#                              #                            legend.position = c(1, .5),
#                              #                            legend.justification = "right",
#                              #                            legend.box.just = "right",
#                              #                            legend.margin = margin(6, 6, 6, 6)
#                              #                            )#+
                                                
#                              #                save_plot(filename=paste("synapseProportions_inactive", sep='_'),plot=barPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)






# # ### HISTOGRAMS ### 






# # Function to create histogram plot
# create_histogram_plot <- function(data,median_data, x_variable, median_var, x_label, filename, lower_bound_x, upper_bound_x) {

#   print(paste0("creating a histogram for: ", x_variable) )
#   bin = (upper_bound_x - lower_bound_x) / 50
#   if(x_variable == "quanta_calc") {quanta_vline = TRUE} else {quanta_vline = FALSE}

#   plot <- ggplot(data, aes(x = get(x_variable), group = Ca_expr, fill = Ca_expr, colour = Ca_expr)) +
#     geom_vline(data = median_data, aes(xintercept = get(median_var), colour = Ca_expr), size = 2, lty = "solid") +
#     geom_freqpoly(aes(y = ..ncount..), binwidth = bin, alpha = 0.7, boundary = 0, size = 1.5) +
#     geom_histogram(aes(y = ..ncount..), binwidth = bin, alpha = 0.4, boundary = 0, position = "identity", colour = "black", size = 0.8) +
#     {if (quanta_vline) geom_vline(xintercept = 1, colour="red",lty="dashed",size=1,alpha=0.9 )}+
#     labs(
#       x = x_label,
#       y = expression(N/N[max]),
#       fill = "",
#       tag = ""
#     ) +
#     {if (color_switch == FALSE) scale_colour_brewer(palette = "Dark2")} +
#     {if (color_switch) scale_colour_manual(labels = new_levels, values = color_override)} +
#     {if (color_switch == FALSE) scale_fill_brewer(palette = "Dark2")} +
#     {if (color_switch) scale_fill_manual(labels = new_levels, values = color_override)} +
#     coord_cartesian(ylim = c(0, 1.1), xlim = c(lower_bound_x, upper_bound_x)) +
#     scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
#     theme_tufte() +
#     facet_grid(Ca_expr ~ ., labeller = label_parsed) +
#     my.theme +
#     theme(
#       legend.position = "none",
#       legend.title = element_blank(),
#       strip.text.y.right = element_blank()
#     )

#   save_plot(
#     filename = paste0(filename, ""),
#     plot = plot,
#     plot_width=plot_height*1.2,plot_height=plot_height*1.5, scale_f = scalefactor,dpi_val = 900)

#   #assign(paste0(tools::to_snake_case(x_variable), "_histo"), plot, envir = .GlobalEnv)
# }

# # Create histogram plots
# create_histogram_plot(light_peak_stat_df, my_medians, "amplitude", "median_amplitude", expression("Peak "*Delta*"F/F"), "histo_amplitude", 0, 10)
# create_histogram_plot(light_peak_stat_df, my_medians, "tau_decay_ms", "median_tau", expression(tau[decay] * " (ms)"), "histo_tau_decay_ms", 0, 200)
# create_histogram_plot(light_peak_stat_df, my_medians, "interSpike_ms", "median_interSpike", expression(Delta * "t (ms)"), "histo_interSpike_ms", 20, 100)
# create_histogram_plot(light_peak_stat_df, my_medians, "t_rise", "median_t_rise", expression(t[rise]*" (ms)"), "histo_interSpike_ms", 0, 40)
# create_histogram_plot(light_peak_stat_df, my_medians, "t_half", "median_t_half", expression(t[1/2] * " (ms)"), "histo_interSpike_ms", 0, 200)
# create_histogram_plot(light_peak_stat_df, my_medians, "t_decay", "median_t_decay", expression(t[decay] * " (ms)"), "histo_interSpike_ms", 0, 400)
# create_histogram_plot(light_peak_stat_df, my_medians, "quanta_calc", "median_quanta", expression(E[Delta*"F/F"] / S[Delta*"F/F"]), "histo_quanta", 0, 20)

        
 



# ### CV plots (VIOLINS) ####


# # Function to create CV plot
# create_CV_plot <- function(data, y_variable, y_label, tag_var, filename,y_upper_bound=1.75) {
#   print(paste0("creating a CV violin plot for: ", y_variable) )

#   x_upper_bound = 4 *1.5
#   jitter_x_range = x_upper_bound/20


#   plot <- ggplot(data, aes(x = groupfact, y = get(y_variable), group = groupfact, colour = Ca)) +
#     geom_sina(aes(group=groupfact, colour=groupfact),size=1.5,alpha=0.5)+
#     geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 0.7,alpha=0.85, colour='black') +
    
#     #geom_sina(alpha = 0.4, size = 1) +
    

#     # stat_summary(geom="errorbar", fun.data=mean_se, width=x_range/40, size=1.5, alpha=1, colour="black")+
#     # stat_summary(geom="line", fun.y=mean, size=1.5, alpha=1, colour="black")+
#     # stat_summary(geom="point", fun.y=mean, size=5, stroke=1.2, shape=21, fill="white",alpha=1, colour="black")+
                            

#     labs(
#       y = y_label,
#       x = expression("["*Ca^'2+'*"] (mM)"),
#       colour = "",
#       tag = tag_var
#     ) +
#     scale_y_continuous(breaks = c(0, 0.5, 1.0,1.5)) + #expand=c(0,0)
#     coord_cartesian(ylim = c(0, y_upper_bound)) +
#     {if (color_switch == FALSE) scale_colour_brewer(palette = "Dark2")} +
#     {if (color_switch) scale_colour_manual(labels = new_levels, values = color_override)} +
#     theme_tufte() +
#     my.theme +
#     theme(legend.position = "none")

#   save_plot(filename = paste0(filename, ""), plot = plot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

#     if(!is.null(tag_var)){assign(paste0(y_variable), plot, envir = .GlobalEnv)}
# }

# # Create CV plots
# create_CV_plot(releaseProb_calc, "CV_amplitude", expression("CV of Peak "*Delta*"F/F"), "B", "CV_amplitude_violin",y_upper_bound = 1.1)
# create_CV_plot(releaseProb_calc, "CV_tau_decay", expression("CV of "*tau[decay]), "C", "CV_tau_violin",y_upper_bound=1.1)
# create_CV_plot(releaseProb_calc, "CV_dt", expression("CV of "*Delta*"t"), "D", "CV_dt_violin",y_upper_bound=1.1)
# create_CV_plot(releaseProb_calc, "CV_t_half", expression("CV of "*t[1/2]), NULL, "CV_t_half_violin")
# create_CV_plot(releaseProb_calc, "CV_t_rise", expression("CV of "*t[rise]), NULL, "CV_t_rise_violin")
# create_CV_plot(releaseProb_calc, "CV_t_decay", expression("CV of "*t[decay]), NULL, "CV_t_decay_violin")



# # ###### average Peaks




#     avgtrace_Plot<-ggplot(clean_peaks_normTime, aes(x=new_normTime, y = dFF,  group=Ca, colour=Ca))+
#                                                     #geom_path(colour="black",size=0.4,alpha=0.2)+
#                                                     stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=2.5,alpha=0.8)+
#                                                     geom_vline(xintercept = 0, lty="longdash",colour='black',size=1)+
                                                    
#                                                     labs( x="Normalized time (s)",
#                                                           y=expression(Delta*"F/F"),
#                                                           #title="Stimulus-evoked iGluSnFR3 activity",
#                                                           colour=expression("[Ca"^{"2+"}*"], mM"),
#                                                                         tag = "A")+
#                                                     {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                                     {if(color_switch)scale_colour_manual(labels=new_levels,values=color_override)}+
#                                                     scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
#                                                     theme_tufte()+
#                                                     my.theme+
#                                                     #facet_grid(~Ca, labeller = label_parsed)+
#                                                     theme(legend.title = element_text(colour="black", size=28, family="sans"),
#                                                         legend.position=c(0.75,0.5))

#                                                     #avg_peaks  <<- avg_PP_tracePlot
#                                                      #       rm(avg_PP_tracePlot
                                                        

#                                                      save_plot(filename=paste0("averagePeaks_basalRP"),plot=avgtrace_Plot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                         

#                      #  averaged_traces <<- avgtrace_Plot
                         

# #### single ROI data


# #### GET ALL THE SPATIAL DATA ##### 



# ###### for loop for plotting spatial data
#         png_dir = "Y:\\Sam/paper1_datasets/singleAP_v2/imageFiles/zstack/zstackOutputs_sliced_SQ"
#         prefix_str = "^MAX_GluSnFR3_SyPhy_"
#         suffix_str = "_0pt5Ca_zstack_lowpwr__FLATTENED.png"
#         all_lut_str = "_0pt5Ca_zstack_lowpwr__FLATTENED_COLOR.png"
#         dir_regex = "GluSnFR3(.*?)region\\d\\d_"
# ### unit test
# ### GluSnFR3_SyPhy_dish06_plate01_region02_0pt5Ca_zstack_lowpwr_
# #subset_df<- PPRplot_data %>% dplyr::filter(vid_key %in% vid_keys_to_save,protocol %in% keep_PP)
# #protocol_vec <- unique(subset_df$protocol[which(subset_df$protocol != "singleAP")])
# #print(length(protocol_vec)) 
# #levs <- c("singleAP","PP60","PP75","PP100","PP150","PP500")
# #ordered_protocols <- protocol_vec[order(match(protocol_vec, levs))]
# #print(ordered_protocols)
# #tag_var = c("F","G","H",
#         #            "I","J","K",
#          #           "L","M","N")
#         #tag_var_subset = tag_vars[6:14]

         
#         #ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d")  

#                     #if(unique(vid_df$vid_key) %in% vid_keys_to_save) { 
#                     #    save_bool = TRUE
#                     #    print(paste0( unique(vid_df$vid_key)," is a region we want to save out. Updated the save_bool value to TRUE" ) )
#                     #} else {
#                     #    save_bool = FALSE
#                     #    #print("This is not a region we want to save out. Updated the save_bool value to FALSE")
#                     #}

#                                         #print(paste0(count, " was value of iteration counter at: Ca = ", unique(Ca_df$Ca), " and protocol = ", unique(protocol_df$protocol) ) )
                                        
#                                         #count = count+1
#                                         #if(count == 3) {
#                                         #guide_bool = TRUE 
#                                         #} else { guide_bool = FALSE}
#                                         #tmp_tag_var = tag_var_subset[count]
#                                         #spat_facet_list
                                                            
#                                         #nam <- paste("spat_facet", count, sep = "")
#                                         #assign(nam, tmp_plot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                
#         print("Taking a crack at generating png plots")
#         subset_df<- releaseProb_calc %>% dplyr::filter(vid_key %in% vid_keys_to_save)

#         ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d") 
#         vid_key_vec <- unique(subset_df$vid_key)
#         vars_to_plot<- c("CV_amplitude",'mean_quanta','releaseProbability_perROI')#, #'release_ratio')#'releaseProbability_perROI')#c("releaseProbability_perROI")
#         save_bool=TRUE
#         guide_bool = TRUE

#         guide_limits = list(c(0,0.5),c(5,20),c(0,1))#,c(0.8,2.0))
#         guide_titles = c(expression(atop(CV[Delta*"F/F"], "1"~Ca^'2+')),expression(atop(E[Delta*"F/F"]/S[Delta*"F/F"], "4"~Ca^'2+')),expression(atop(P[iGlu], "0.5"~Ca^'2+')))#, expression(frac(E/S["4Ca"], E/S["2Ca"])))    #expression(P[iGlu], "0.5 mM"~Ca^'2+'))# # expression(CV[tau*"decay"])#
        
#         #expression("Avg."~E[Delta*"F/F"] / S[Delta*"F/F"])
#         tmp_tag_var_override = c("A","B","C")#c("E","F", "G")


#         ###scale bar params! set the left-corner 
#         user_x = 20
#         user_y = 24.86




#         for (k in 1:length(vid_key_vec)){
            
#                 vid_df<- releaseProb_calc %>% dplyr::filter(vid_key == vid_key_vec[k]) %>% mutate(fix_Ca = paste0(Ca_mM, "Ca") )
#                 Ca_vec <- unique(vid_df$fix_Ca)
#                 print(Ca_vec)
        
#                 max_gradient_val = max(vid_df$releaseProbability_perROI)

                
#                     # if(!is.null(vid_df) & save_bool==TRUE){
                    
#                     #    png_plotter_v4(df = vid_df, png_dir = png_dir,prefix_str=prefix_str, suffix_str = suffix_str,
#                     #                     all_lut_str = all_lut_str, dir_regex = dir_regex, vars_to_plot=vars_to_plot,save_bool=save_bool,
#                     #                     tag_var = tmp_tag_var_override[1],user_x = user_x, user_y = user_y)

#                     #    png_plotter_ROIoverlay(df=vid_df,png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,
#                     #                             dir_regex=dir_regex,vars_to_plot=vars_to_plot,save_bool=save_bool,
#                     #                             tag_var = tmp_tag_var_override[2],ROI_IDs = ROInames,user_x = user_x, user_y = user_y)
                                
#                     # }     

#                 count=0
#                # for(i in 1:length(Ca_vec)) {
#                                     #tmp_tag_var = "I"
#                                     tmp_plot = NULL
#                                     #print("printing 0.5 mM Ca df")
                                    

#                                     #Ca_df_0pt5<- vid_df %>% dplyr::filter(fix_Ca == "0.5Ca")
#                                     #Ca_df_0pt5$releaseProbability_perROI[Ca_df_0pt5$releaseProbability_perROI == 0] <- NA
                                    

#                                     #print(Ca_df_0pt5)
#                                     Ca_df_1<- vid_df %>% dplyr::filter(fix_Ca == "1Ca")
#                                     Ca_df_4<- vid_df %>% dplyr::filter(fix_Ca == "4Ca")
#                                     Ca_df_05<- vid_df %>% dplyr::filter(fix_Ca == "0.5Ca")
#                                     #Ca_df_saturated <- vid_df %>% dplyr::filter(fix_Ca == "4Ca")
#                                     #Ca_df_saturated<- left_join(Ca_df_saturated, corr_stats_saturation_score)
#                                     Ca_df_list = list(Ca_df_1,Ca_df_4,Ca_df_05)#, Ca_df_saturated)

#                                     for(j in 1:length(vars_to_plot)){

#                                             #check_data<- Ca_df %>% select(fix_Ca, vid_key,releaseProbability_perROI)
#                                             #print("This is the current dataset:")
#                                             #print(check_data)
#                                         if(save_bool==TRUE) { 
#                                             guide_bool = TRUE
#                                             var_to_plot = vars_to_plot[j]
#                                             tmp_guide_title = guide_titles[j]
#                                             tmp_limits = guide_limits[[j]]
#                                             tmp_tag = tmp_tag_var_override[j]
#                                             if(j == 1){tmp_scalebar = TRUE} else {tmp_scalebar = FALSE}
#                                             Ca_df = Ca_df_list[[j]]
#                                             print(Ca_df)

#                                            tmp_plot = png_plotter(df=Ca_df,
#                                                                     png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
#                                                                     vars_to_plot=var_to_plot, 
#                                                                     save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
#                                                                     ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar)
#                                     }

#                                     }

                                                                    

#                 #}

#         }
                                    
                                        



# ##### step 1. establish subset of the dataframe that can be used to plot traces and physiology

# double_check<- clean_df %>% dplyr::filter(trackROI_key %in% ROIs_to_save) %>% ungroup() %>% group_by(ROI_key) %>% 
#                         dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5) ) %>%
#                         summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
#                         dplyr::filter(maxInterframe > 0.2) %>%
#                         select(ROI_key,maxInterframe) 

# print(double_check)
#                         #dplyr::filter(maxInterframe > 0.1)
# print(paste0("These are the bad recordings: ", unique(double_check$ROI_key)))


# tmp_trace_data <- clean_df %>% dplyr::filter(trackROI_key %in% ROIs_to_save)#, !ROI_key %in% unique(double_check$ROI_key))
# print(unique(tmp_trace_data))
# tmp_timezero_df<- suppressMessages(tmp_trace_data %>%  group_by_at(groupers) %>% summarise(time_zero = unique(firstStim[!is.na(firstStim)] ) ) )
# tmp_trace_data<- suppressMessages(left_join(tmp_trace_data,tmp_timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))

# tracked_ROIs <- unique(tmp_trace_data$trackROI_key)
# tmp_tag_vars <- list(c("D","G"),#c("H","K"), #trace
#                     c("E","H"),#c("I", "L"), #releaseprob
#                     c("F","I"))#c("J","M"))#, #CV amplitude
#                     #c("K", "N")) #amplitude
# for (i in 1:length(tracked_ROIs)) {
#                     count = i
#                     print(paste0("loop iter on count: ", i))
#                     #print("Generating data visualizations for trackROI_key:", ROI)
#                     ROI_title<-gsub("GluSnFR3-dish\\d\\d-plate\\d\\d-region\\d\\d-","",tracked_ROIs[i])


#                     tmp_original_stats<- peak_stat_df %>% 
#                                          dplyr::filter(trackROI_key == tracked_ROIs[i])
#                     tmp_summary_data <- releaseProb_calc %>%
#                                         dplyr::filter(trackROI_key == tracked_ROIs[i])
#                     tmp_ROIs_with_peaks<- ROIs_with_peaks %>% 
#                                         dplyr::filter(trackROI_key == tracked_ROIs[i],
#                                                         peak_positive > 0)
#                     positive_ROI_list<- unique(tmp_ROIs_with_peaks$ROI_key)


#                     print(tmp_summary_data)

#                     tmp_trace_subset<- tmp_trace_data %>% 
#                                         dplyr::filter(trackROI_key == tracked_ROIs[i])#,
#                                                         #ROI_key %in% positive_ROI_list)
#                     tmp_trace_subset$Ca<- as.factor(as.character( (tmp_trace_subset$Ca) ) )
#                     levels(tmp_trace_subset$Ca) <- new_levels


                        
#                     derived_amplitudes<- left_join(tmp_trace_subset,light_spont_avg) %>% 
#                                          group_by(Ca, Ca_mM, trackROI_key,ROI_key) %>%
#                                          dplyr::filter(new_normTime > 0, new_normTime < 0.5, peakID != "NotPeak") %>%
#                                          summarise(amplitude = max(dFF, na.rm=TRUE),
#                                                     quanta_calc = amplitude/spont_amplitude) %>%
#                                          dplyr::filter(trackROI_key == tracked_ROIs[i])
#                     derived_amplitudes$groupfact <- factor(derived_amplitudes$Ca_mM)    

#                     max_amp = derived_amplitudes %>% dplyr::filter(Ca_mM == 4)
#                     max_amp = max(max_amp$amplitude)
                    
#                     print(paste0("The maximum amplitude for :",unique(derived_amplitudes$trackROI_key), "was ", max_amp," DeltaF/F"))


#                     tmp_amplitude<-ggplot(derived_amplitudes, aes(x=Ca_mM, y=quanta_calc))+#linetype=segmentation,
#                                              #geom_sina(aes(group=groupfact, colour=groupfact),size=2,alpha=0.5)+
                            
#                                              geom_point(aes(group=groupfact,colour=groupfact),alpha=0.5, size=2)+
                                             
#                                              geom_hline(yintercept = 1, size=0.7,lty="dashed",colour="black")+
#                                              stat_summary(geom="line", fun.y = mean, size=1, alpha=0.85,colour="black" )+
                                             
#                                              stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=0.8, alpha=1, colour="black")+
#                                              stat_summary(geom="point", fun.y=mean, color="black",shape=21,stroke=1.1, size=2.5,fill="white",alpha=1)+
                                             
#                                              labs( y=expression(E[Delta*"F/F"] / S[Delta*"F/F"]),
#                                                     x=expression("["*Ca^'2+'*"] (mM)"),
#                                                     title="",
#                                                     colour="",
#                                                     tag = tmp_tag_vars[[3]][i]
#                                                     )+
#                                             {if(color_switch)scale_colour_manual(values=c("0.5" = color_override[1],"1" = color_override[2],"2" = color_override[3],"4" = color_override[4]) )}+
#                                             #{if(color_switch)scale_alpha_manual(values=c(0.5,0.5,0.5,0.5))}+
                                                    
#                                             scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
#                                             scale_y_continuous(breaks=c(1,10,20,30))+
#                                             coord_cartesian(ylim=c(0,35),xlim=c(0.25,4.25))+
                            
#                                             theme_tufte()+
#                                             my.theme+
#                                             theme(legend.position = "none",
#                                                   #axis.text.x=element_text(colour="black", size=28*scalefactor, family="sans",angle=45, hjust=1)
#                                                   )
                    
#                     nam <- paste("tmp_amplitude", count, sep = "")
#                     assign(nam, tmp_amplitude,envir = .GlobalEnv)
                    
                      

#                     tmp_CV_amplitude<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_amplitude))+#linetype=segmentation,
#                                               geom_line(size=1,color="black",alpha=0.85)+
#                                                 geom_point(color="black",shape=21,stroke=1.1, size=2.5,fill="white")+
                                                                 
#                                                 labs( x=expression("["*Ca^'2+'*"] (mM)"),
#                                                         y=expression("CV"[Delta*F/F]),
#                                                         title="",
#                                                         tag = tmp_tag_vars[[3]][i]
#                                                         )+
#                                                 {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                                 {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
#                                                 #guides(colour="none")+
#                                                 scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
#                                                 scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
#                                                 coord_cartesian(ylim=c(0,0.5),xlim=c(0.25,4.25))+
#                                                                #, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
#                                             theme_tufte()+
#                                             my.theme+
#                                             theme(legend.position = "none")

#                         nam <- paste("tmp_CV_amplitude", count, sep = "")
#                         assign(nam, tmp_CV_amplitude,envir = .GlobalEnv)
                                
                        


#                         tmp_CV_tau_decay<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_tau_decay))+#linetype=segmentation,
#                                              geom_line(size=2,color="black")+
#                                              geom_point(size=4,color="black")+
#                                              labs( y=expression("CV of "*tau[decay]*", ms"),
#                                                     x=expression("[Ca"^{"2+"}*"], mM"),
#                                                     title="",
#                                                     colour=""
#                                                     )+
#                                             coord_cartesian(ylim=c(0,1.0))+
#                                             scale_x_continuous(breaks=c(0.5,1,2,4))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
#                                             theme_tufte()+
#                                             my.theme+
#                                             theme(legend.position = "none")
                                            
                       



#                         tmp_CV_dt<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=CV_dt))+#linetype=segmentation,
#                                              geom_line(size=2,color="black")+
#                                              geom_point(size=4,color="black")+
#                                              labs( y=expression("CV of "*Delta*"t, ms"),
#                                                     x=expression("[Ca"^{"2+"}*"], mM"),
#                                                     title="",
#                                                     colour=""
#                                                     )+
#                                             coord_cartesian(ylim=c(0,1.0))+
#                                             scale_x_continuous(breaks=c(0.5,1,2,4))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
#                                             theme_tufte()+
#                                             my.theme+
#                                             theme(legend.position = "none")
                                            
                      

#                         # tmp_N_fit<- tmp_summary_data #%>% dplyr::filter(Ca_mM == 4)
#                         # y_limit = ifelse(max(unique(tmp_summary_data$var_amplitude),na.rm=TRUE) > 1, max(unique(tmp_summary_data$var_amplitude),na.rm=TRUE), 1)  
#                         # add_text_labels = TRUE
#                         # #y_limit = 5
#                         # #### try fitting an individual point? 
#                         # tmp_N_sites<-ggscatter(tmp_N_fit, 
#                         #                     x = "mean_amplitude", y = "var_amplitude", 
#                         #                     size = 3,
#                         #                     color='Ca',
#                         #                     shape=21,
#                         #                     stroke=1.1,
#                         #                     fill='white',
#                         #                     alpha=1#,
#                         #                     )+
#                         #     #geom_vline(xintercept = 1, lty='dashed',size=0.5,colour="black")+
                            
#                         #     {if(add_text_labels)stat_regline_equation(label.x = 1, label.y=0.95*y_limit,size=7,formula = y ~ poly(x, 2),
#                         #                             aes(label =  paste(..eq.label.., sep = "~~~~")))} +
#                         #     {if(add_text_labels)geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)}+
#                         #     #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..r.label.., sep = "~")), label.x = 1, label.y=0.95*max_y, size=7, geom = "text")} +
#                         #     #{if(add_text_labels)ggpubr::stat_cor(aes(label = paste(..p.label.., sep = "~")), label.x = 1, label.y=0.82*max_y, size=7, geom = "text")} +
                            
#                         #     #{if(add_text_labels)stat_regline_equation(label.x = 1, label.y=0.82*max_y,size=7)}+
                              
#                         #     #scale_fill_manual(values = c("red"  = "white",'blue' = 'blue'), labels = c("red" = "More Glutamate Released", "blue" = "Saturated"))+
#                         #     scale_colour_manual(values = color_override)+     
#                         #     coord_cartesian(xlim=c(0,10),ylim=c(0,y_limit))+
#                         #     #scale_y_continuous(breaks=c(0,0.5,1.0))+

#                         #     #coord_flip()+
#                         #     labs(x=expression(Delta*"F/F"),
#                         #             y=expression(sigma[Delta*"F/F"]^2),
#                         #             tag=""
#                         #             )+
#                         #     theme_tufte()+
#                         #     #facet_grid(~RRP_score, labeller = labeller(RRP_score = RRP_score.labs))+
#                         #     my.theme+
#                         #     theme(legend.position = "none")
                            
#                         #     nam <- paste("tmp_N_sites", count, sep = "")
#                         #     assign(nam, tmp_N_sites,envir = .GlobalEnv)
                    


#                         tmp_releaseProb<- ggplot(tmp_summary_data, aes(x=Ca_mM, y=releaseProbability_perROI))+
                                                
#                                                 geom_line(size=1,color="black",alpha=0.85)+
#                                                 geom_point(color="black",shape=21,stroke=1.1, size=2.5,fill="white")+
                                                                 
#                                                 labs( x=expression("["*Ca^'2+'*"] (mM)"),
#                                                         y=expression(P[iGlu]),
#                                                         title="",
#                                                         tag = tmp_tag_vars[[2]][i]
#                                                         )+
#                                                 {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                                 {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
#                                                 #guides(colour="none")+
#                                                 scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
#                                                 scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
#                                                 coord_cartesian(ylim=c(0,1.0),xlim=c(0.25,4.25))+
                                                                    
#                                                 theme_tufte()+
#                                                 my.theme+
#                                                 theme(legend.position = "none",
#                                                         #axis.text.x=element_text(colour="black", size=28*scalefactor, family="sans",angle=45, hjust=1)
#                                                         )
#                         nam <- paste("tmp_releaseProb", count, sep = "")
#                         assign(nam, tmp_releaseProb,envir = .GlobalEnv)
                    


#                         tmp_display_traces<- tmp_trace_subset #%>% dplyr::filter(!ROI_key %in% unique(double_check$ROI_key))


                        # tmp_tracePlot<- ggplot(tmp_display_traces, aes(x=new_normTime, y = dFF,  group=1, colour=Ca,alpha=Ca))+
                        #                             geom_vline(xintercept = 0, lty="longdash",colour='grey21',size=0.9)+
                                                    
                        #                             geom_path(size=0.5)+
                        #                             stat_summary_bin(aes(group=Ca), colour="black",geom="line", fun=mean, binwidth=interFrame.var,size=1,alpha=0.95)+
                                                    
                        #                             labs( x="Normalized time (s)",
                        #                                   y=expression(Delta*"F/F"),
                        #                                   title=paste0("Bouton",gsub("ROI00", "  ", ROI_title)),
                        #                                   #colour=expression("[Ca"^{"2+"}*"], mM"),
                        #                                                 tag = tmp_tag_vars[[1]][i])+
                        #                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                        #                             {if(color_switch)scale_colour_manual(labels=new_levels,values=color_override)}+
                        #                             {if(color_switch)scale_alpha_manual(labels=new_levels,values=c(0.8,0.7,0.5,0.35))}+
                        #                             coord_cartesian(ylim=c(0,12.7))+
                        #                             scale_y_continuous(breaks=c(0,2.5,5,7.5,10,12.5))+
                        #                             scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
                        #                             theme_tufte()+
                        #                             my.theme+
                        #                             facet_grid(~Ca, labeller = label_parsed)+
                        #                             theme(legend.title = element_text(colour="black", size=28, family="sans"),
                        #                                 legend.position="none",
                        #                                 strip.text.x = element_text(colour="black",size=30*scalefactor,family="sans"))

                                                        

#                                                      ggsave(filename=paste0("avg_traces_", tracked_ROIs[i],"_"),plot=tmp_tracePlot, device="jpeg",dpi=600, bg="white",units="in",width=plot_height*1.8,height=plot_height)
                                                         
#                                                        if(count <= 2){
#                                                         nam <- paste("tmp_tracePlot", count, sep = "")
#                                                         assign(nam, tmp_tracePlot,envir = .GlobalEnv)
#                                                         }

#                         inset_traces = tmp_display_traces %>% dplyr::filter(Ca_mM == 0.5)
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
                      


#                         # gs = list(tmp_tracePlot,  #1
#                         #          tmp_releaseProb,
#                         #          tmp_N_sites,       #2
#                         #          tmp_amplitude)  #3
#                         # margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
                    



#                         # hlay <- rbind(c(1,1,1,1,1,1,2,2,2,3,3,3),
#                         #                c(1,1,1,1,1,1,2,2,2,3,3,3),
#                         #                c(1,1,1,1,1,1,2,2,2,3,3,3),
                                       
#                         #                c(1,1,1,1,1,1,4,4,4,4,4,4),
#                         #                c(1,1,1,1,1,1,4,4,4,4,4,4),
#                         #                c(1,1,1,1,1,1,4,4,4,4,4,4)
#                         #            )
                
                    

                    
#                         # ROI_report = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
#                         # save_plot(filename=paste0("ROI_report_",tracked_ROIs[i],"_"),plot=ROI_report, plot_width=28,plot_height=12, scale_f = scalefactor,dpi_val = 600)



# }

                               

        





# list.output = list(ROI_totals = getROIs,tracked_ROI_totals = get_trackedROIs, ROI_probability = ROIs_with_peaks, releaseProb_calc = releaseProb_calc, Medians = my_medians,Averages = plot_avgs) #, Averages = plot_avgs, Original_stats =peak_stat_df)
# list.output
#full_model_output = N_sites_reg_model, reasonable_N_sites = check_good_N_sites, fit_predict = generate_fit_predict, sum_sq_resid = get_sum_sq_resid, 
list.output = list(clean_df = clean_df, ROIs_with_peaks = ROIs_with_peaks, light_peak_stat_df=light_peak_stat_df, join_df = join_df, releaseProb_calc = releaseProb_calc,corr_stats = corr_stats,data_for_umap =data_for_umap)#umap_data = df_annotated, scaled_df = standardized_umap) #ROIs_with_peaks = ROIs_with_peaks, silentROIs = silent_ROIs, true_noise = true_noise_amplitude)#
}
						




