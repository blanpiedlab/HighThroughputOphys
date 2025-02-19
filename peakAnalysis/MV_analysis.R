
subDir <- "MV_analysis_v2"   


source(paste0(path,"sub-functions/setFigureDirectory.R") ) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"sub-functions/def_stim_vlines.R"))
source(paste0(path,"peakAnalysis/find_replicates.R"))
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))




library(ggforce)
library(scales)
library(nls.multstart)

library(umap)




MV_analysis<- function(peak_stat_df = peak_stats, 
                       df = df_traces, 
                       spont_avgs = spont_avg_df, 
                       groupers = groupers, 
                       color_override = color_override, 
                       remove_cols=drops, 
                       plotBy = plotBy, 
                       levels=levels,
                       remove_col = drops,
                       ROIs_to_save = ROIs_to_save){



############################## figure out new mean_variance analysis. Copy pasted wholesale from releaseProb_v2.R to attempt to make more modular code. ########################################

                            
                              #could be deprecated
                              vid_keys_to_save = unique(str_extract(ROIs_to_save, "dish\\d\\d-plate\\d\\d-region\\d\\d") ) 
                              
                              #essential parameter definitions
                              color_switch = !is.null(color_override)
                              plot_height = 8
                              interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE) 


                              add_second_model = FALSE

                              

                              #essential parameter
                              remove_neurons = c("dish05-plate03-region02", "dish05-plate06-region01", "dish05-plate06-region02")#,"dish05-plate01-region01")
                              

                              #### clean up dataframe by merging with peakStats. could be redundant

                              peak_stat_df <- peak_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID), 
                                                                      Ca_mM = case_when(Ca == "0pt5Ca" ~ 0.5,
                                                                                        Ca == "1Ca" ~ 1,
                                                                                        Ca == "2Ca" ~ 2,
                                                                                        Ca == "4Ca" ~ 4),
                                                                      imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>%
                                                        dplyr::filter(!imaging_region %in% remove_neurons)
                              
                              light_df <- df[,!(names(df) %in% remove_col)]



                              #essential to eliminate cross-talk with bad recordings
                              ROIs_to_remove <- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are finite
                                                       group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%
                                                       dplyr::filter( is.na(dFF))
                              ROIs_to_remove_timeskip<- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are do not have a frame skip in time of interest
                                                       group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%            
                                                       summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
                                                       dplyr::filter(maxInterframe > 0.18) %>% #empirically defined value which eliminates all bad recordings for subsequent visualizations
                                                       select(ROI_key,maxInterframe) 



#uncomment these lines to get a readout of ROI removal step in the terminal                           
#print(paste0("These are the bad recordings: ", unique(ROIs_to_remove$ROI_key)))
#print(paste0("These are the bad recordings with a frame jump during post-stimulus: ", unique(ROIs_to_remove_timeskip$ROI_key)))

                       
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

                        #might be a premature removal of these data
                        #rm(clean_peaks,timezero_df)
                        
                        clean_peaks_normTime<- suppressMessages(left_join(clean_peaks,timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))
                        

                        #refactor the dataframe for appropriate facet_labels
                        new_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
                        
                        clean_peaks_normTime$Ca<- as.factor(as.character( (clean_peaks_normTime$Ca) ) )
                        peak_stat_df$Ca<- as.factor(as.character( (peak_stat_df$Ca) ) )

                        levels(clean_peaks_normTime$Ca) <- new_levels
                        levels(peak_stat_df$Ca) <- new_levels
                        
                        
                        
                      

                        #### generate ROI count from raw data, could be deprecated

                        getROIs<- clean_df %>% group_by(vid_key,Ca, Ca_mM) %>% summarise(ROI_count = length(unique(ROI_key)))
                        get_trackedROIs <- clean_df %>% group_by(Ca) %>% mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key) ) %>% 
                                                                        summarise(imaging_region_count = length(unique(imaging_region)),
                                                                                    trackedROI_count  = length(unique(trackROI_key)))                        

                        area_imaging_region_um = 25.6 * 25.6

                        #### Statistics about imaging regions for visualization, could be deprecated
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






                    # Defining Could be deprecated
                    ROIs_with_peaks<- suppressMessages(clean_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5), !ROI_key %in% ROIs_to_remove) %>%  
                                                                                group_by(Ca, Ca_mM, vid_key,trackROI_key,ROI_key,ROINumber,exposeNum,replicates) %>%  
                                                                                summarise(peak_detected = length(unique(peakID))-1,
                                                                                                peak_positive = case_when(peak_detected > 0 ~ 1,
                                                                                                                          peak_detected == 0 ~ 0) ) 
                                                                                )

                    droppers = c("sensor","marker","segmentation","TTL_start","dish","plate","region","peakID","protocol",'exposeNum', "windowedPeakID","Ca")
                    cleanup_peak_stat_data_point <- peak_stat_df %>% dplyr::filter(ROI_key == "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005", interSpike_ms < 100)
                    cleanup_peak_stat_df <- peak_stat_df %>% dplyr::filter(ROI_key != "GluSnFR3-SyPhy-marker-dish05-plate01-region04-0pt5Ca-singleAP-repl07-ROI0005" )
                    cleanup_merge <- bind_rows(cleanup_peak_stat_df, cleanup_peak_stat_data_point)

                    light_peak_stat_df <- cleanup_merge[,!(names(cleanup_merge) %in% droppers)]
                    rm(cleanup_peak_stat_data_point, cleanup_peak_stat_df, cleanup_merge)
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
                    
                    ## refactor
                    light_peak_stat_df$Ca_expr <- factor(light_peak_stat_df$Ca_expr, levels = c("0pt5Ca","1Ca","2Ca","4Ca"))
                    levels(light_peak_stat_df$Ca_expr) <- new_levels

                                                            
                    

                    




                    estimate_max_signal_df<- clean_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.15)) %>%
                                                         group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>% 
                                                         summarise(amplitude = max(dFF, na.rm=TRUE))

                    mean_noise_df <- clean_df %>% dplyr::filter(peakID == "NotPeak") %>%
                                                         group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>% 
                                                         summarise(median_noise = median(dFF, na.rm=TRUE))

                    estimates_per_trial_df <- left_join(estimate_max_signal_df,mean_noise_df) %>% mutate(normalize_amplitude = amplitude - median_noise)
                    estimates_per_trial_df <- left_join(estimates_per_trial_df, light_spont_avg) %>% mutate(corrected_amplitude = case_when(normalize_amplitude < spont_amplitude ~ 0,
                                                                                                                                             normalize_amplitude >= spont_amplitude ~ normalize_amplitude),
                                                                                                            quanta_calc = corrected_amplitude/spont_amplitude)

                    MV_per_synapse_df <- estimates_per_trial_df %>% group_by(Ca, Ca_mM, trackROI_key) %>% 
                                                                    summarise(n_obs = n(),
                                                                              mean_amplitude = mean(corrected_amplitude,na.rm=TRUE),
                                                                              var_amplitude = var(corrected_amplitude,na.rm=TRUE),
                                                                              mean_quanta = mean(quanta_calc, na.rm=TRUE), 
                                                                              var_quanta = var(quanta_calc,na.rm=TRUE),
                                                                              sd_amplitude = sd(corrected_amplitude,na.rm=TRUE),
                                                                              se_amplitude = sd_amplitude/sqrt(n_obs),
                                                                              CV_amplitude = sd_amplitude/mean_amplitude

                                                                              ) 
                    get_quanta_scores_4Ca<- MV_per_synapse_df %>% dplyr::filter( Ca == "4Ca") %>% ungroup() %>% select(trackROI_key, mean_quanta, var_quanta)
                    get_quanta_scores_2Ca<- MV_per_synapse_df %>% dplyr::filter( Ca == "2Ca") %>% ungroup() %>% select(trackROI_key, mean_quanta, var_quanta)

                    names(get_quanta_scores_4Ca)[names(get_quanta_scores_4Ca) == 'mean_quanta'] <- 'mean_quanta_4Ca'
                    names(get_quanta_scores_4Ca)[names(get_quanta_scores_4Ca) == 'var_quanta'] <- 'var_quanta_4Ca'
                    names(get_quanta_scores_2Ca)[names(get_quanta_scores_2Ca) == 'mean_quanta'] <- 'mean_quanta_2Ca'
                    names(get_quanta_scores_2Ca)[names(get_quanta_scores_2Ca) == 'var_quanta'] <- 'var_quanta_2Ca'
                    
                    MV_per_synapse_df <- left_join(MV_per_synapse_df,get_quanta_scores_2Ca)
                    MV_per_synapse_df <- left_join(MV_per_synapse_df, get_quanta_scores_4Ca) %>% mutate(quanta_ratio = mean_quanta_4Ca/mean_quanta_2Ca,
                                                                                                         var_ratio = var_quanta_4Ca/var_quanta_2Ca,
                                                                                                         MV_score = case_when(quanta_ratio > 1 & var_ratio > 2 ~ "crosstalk?",
                                                                                                                              quanta_ratio < 1 & var_ratio <= 1 ~ "exhausted terminal?",
                                                                                                                              quanta_ratio > 1 & var_ratio <= 2 ~ "MV behavior"))
                    print(MV_per_synapse_df) 
                                                                     
                    keep_ROIs_MV <- MV_per_synapse_df %>% dplyr::filter( Ca == "4Ca")
                    keep_ROIs_MV_list<- unique(keep_ROIs_MV$trackROI_key)
                    MV_per_synapse_df <- MV_per_synapse_df %>% dplyr::filter(trackROI_key %in% keep_ROIs_MV_list)




                    #uncomment to only sample a subset of fits
                    # n=50 
                    # sample_ROIs<- unique(MV_per_synapse_df$trackROI_key)
                    # randomROIs<- sample(sample_ROIs,n)
                    # MV_per_synapse_df <- MV_per_synapse_df %>% dplyr::filter(trackROI_key %in% sample_ROIs) #%>% dplyr::filter(Ca != "4Ca")
                     


                    # ####generate releaseProb_calc version of N_sites

                    library(nls.multstart) 



                    N_sites_reg_model<- MV_per_synapse_df %>% ungroup() %>% group_by(trackROI_key) %>% 
                                            select(trackROI_key, mean_quanta,var_quanta) %>% 
                                            #dplyr::filter(!is.na(var_amplitude), !is.na(mean_amplitude), trackROI_key %in% active_4Ca_list, trackROI_key %in% active_2Ca_list, trackROI_key %in% active_1Ca_list) %>%
                                            nest() %>%
                                             mutate(binom_fit = purrr::map(data, ~nls_multstart(var_quanta ~ Q*mean_quanta - mean_quanta^2/N,
                                                                            data=.x,
                                                                            iter =1000,
                                                                            start_lower = c(Q = 0.7, N = 1),
                                                                            start_upper = c(Q = 1.3, N = 15),
                                                                            supp_errors = 'Y',
                                                                            na.action = na.omit,
                                                                            lower = c(Q = 0.7, N= 1),
                                                                            upper = c(Q = 1.3, N = 100))),
                                                    tidied = map(binom_fit, tidy),
                                                    augment = map(binom_fit, augment)
                                                    ) 

                    check_good_N_sites <- N_sites_reg_model %>% unnest(tidied) %>%
                                            select(trackROI_key,term,estimate) %>%
                                            spread(term,estimate) %>%
                                            mutate(N_min = N) %>%
                                        mutate(N_sites_group = ifelse(N_min < 20, "A", "B")) %>%
                                        dplyr::filter(N_min > 0, N_min < 60)

                    print(check_good_N_sites[,c(1:5)])

                    good_N_sites<- MV_per_synapse_df 
                    good_N_sites <- left_join(check_good_N_sites, good_N_sites) %>% mutate(P_rw = mean_quanta/(N_min*Q),
                                                                                            Q_w = Q)


                    generate_fit_predict<- check_good_N_sites %>% group_by(across(everything())) %>%
                                                                  summarise(amp_points = seq(from = 0, to = 30, length.out=300)) %>%
                                                                  ungroup() %>%
                                                                  group_by(trackROI_key) %>% 
                                                                  mutate(var_pred = Q*amp_points - amp_points^2/N,
                                                                         N_sites_group = ifelse(N_min < 20, "A", "B"))




                    generate_fit_predict <- generate_fit_predict %>% group_by(trackROI_key) %>% select(trackROI_key,  amp_points, var_pred) #N_sites_group,
                    names(generate_fit_predict)[names(generate_fit_predict) == 'amp_points'] <- 'mean_quanta'
                    names(generate_fit_predict)[names(generate_fit_predict) == 'var_pred'] <- 'var_quanta'

                    predict_apex<- generate_fit_predict %>% group_by(trackROI_key) %>% summarise(max_quanta_fit = mean_quanta[which.max(var_quanta)])
                    compare_fit_apex <- MV_per_synapse_df %>% ungroup() %>% group_by(trackROI_key) %>% summarise(max_quanta = max(mean_quanta, na.rm=TRUE)) 
                    compare_fit_apex<- left_join(compare_fit_apex, predict_apex) %>% group_by(trackROI_key) %>% mutate(binom_fit_score = ifelse(max_quanta > max_quanta_fit,TRUE,FALSE))
                    print(compare_fit_apex) 
                    fittable_ROIs<- compare_fit_apex %>% dplyr::filter(binom_fit_score == TRUE)


                    check_fitted<-N_sites_reg_model  %>% unnest(augment)
                    get_sum_sq_resid <- check_fitted %>% group_by(trackROI_key) %>%
                                                      summarise( sum_sq_resid = sum(.resid^2,na.rm=TRUE))

                    keep_ROIs  <- get_sum_sq_resid %>% dplyr::filter(sum_sq_resid <= 1)  


if(add_second_model == TRUE) {

                    N_sites_nonunif_reg_model<- MV_per_synapse_df %>% ungroup() %>% group_by(trackROI_key)  %>%
                                            select(trackROI_key, mean_quanta,var_quanta) %>% 
                                            nest() %>%
                                             mutate(binom_fit = purrr::map(data, ~nls_multstart(var_quanta ~ ( ( Q*mean_quanta - (Q*mean_quanta^2 * (1 + alpha) )/(mean_quanta + N*Q*alpha) )) ,
                                                                            data=.x,
                                                                            iter = 1000,
                                                                            start_lower = c(Q = 0.7, N = 1,alpha=0),
                                                                            start_upper = c(Q = 1.3, N = 15,alpha=1000),
                                                                            supp_errors = 'Y',
                                                                            na.action = na.omit,
                                                                            lower = c(Q = 0.7, N= 1, alpha=-Inf),
                                                                            upper = c(Q = 1.3, N = 100, alpha=Inf))),
                                                    tidied = map(binom_fit, tidy),
                                                    augment = map(binom_fit, augment)
                                                    ) #%>%
                                            #unnest(tidied) %>%
                                            #select(trackROI_key,term,estimate) %>%
                                            #spread(term,estimate) %>%

                    check_nonunif_model <- N_sites_nonunif_reg_model %>% unnest(tidied) %>%
                                            select(trackROI_key,term,estimate) %>%
                                            spread(term,estimate) #%>%
                                            #mutate(Pr = mean_amplitude/(N*Q))

                    nonunif_N_sites<- MV_per_synapse_df
                    nonunif_N_sites<- left_join(check_nonunif_model,nonunif_N_sites) #%>% mutate(Pr = mean_amplitude/(N*Q))

                    generate_nonunif_fit_predict<- check_nonunif_model %>% group_by(across(everything())) %>%
                                                                  summarise(amp_points = seq(from = 0, to = 30, length.out=300)) %>%
                                                                  ungroup() %>%
                                                                  group_by(trackROI_key) %>% 
                                                                  mutate(var_pred = ( ( Q*amp_points - (Q*amp_points^2 * (1 + alpha) )/(amp_points + N*Q*alpha) )) ) 

                    names(generate_nonunif_fit_predict)[names(generate_nonunif_fit_predict) == 'amp_points'] <- 'mean_quanta'
                    names(generate_nonunif_fit_predict)[names(generate_nonunif_fit_predict) == 'var_pred'] <- 'var_quanta'

                    nonunif_N_sites <- nonunif_N_sites %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key))  #%>% dplyr::filter(N_min >= 0 )
                    generate_nonunif_fit_predict <- generate_nonunif_fit_predict %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key))
                    get_sum_sq_resid_nonunif <- N_sites_nonunif_reg_model %>% unnest(augment) %>% group_by(trackROI_key) %>% summarise( sum_sq_resid = sum(.resid^2,na.rm=TRUE))

                    nonunif_N_sites<- nonunif_N_sites %>% mutate(groupfact = "nonUnif")
                    releaseProbPlot<-ggplot(nonunif_N_sites, aes(x=groupfact,y=N))+
                            geom_sina(size=2,alpha=0.5, colour="magenta")+
                            geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 1.2,alpha=0.85, colour='black') +
    
                            labs( x=expression("placeholder"),
                                    y=expression(N[sites]),
                                    #title=expression("Binom. Uniform "*P[r]),
                                    tag= "E"#"F" Figure 3F in manuscript, changed for SFN
                                    )+
                            
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")

                     save_plot(filename = "nonunif_model_N", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                     #Pr_distribution<<- releaseProbPlot


}                 

                   # model_outputs<- full_join(good_N_sites,nonunif_N_sites)  
good_N_sites<- good_N_sites %>% mutate(model_group = "Unif")


  disp_N_sites<- good_N_sites %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key)) 


  Pr_filter<- disp_N_sites %>% ungroup() %>% group_by(trackROI_key) %>% summarise(max_Pr = max(P_rw, na.rm=TRUE)) %>%
                                             ungroup() %>%
                                             dplyr::filter(max_Pr > 0.45)

  disp_N_sites <- disp_N_sites %>% dplyr::filter(trackROI_key %in% unique(Pr_filter$trackROI_key))
 

  releaseProbPlot<-ggplot(disp_N_sites, aes(x=model_group, y=N))+
                            geom_sina(size=2,alpha=0.5, colour="grey41")+
                            geom_violin(fill = NA, draw_quantiles = c(0.25, 0.5, 0.75), size = 1.1,alpha=0.85, colour='black') +
    
                            labs( x=expression(N[sites]),
                                    y=expression(N[sites]),
                                    subtitle=str_wrap("Estimated SVs released per stimulus",25),
                                    tag="E" #"G"Figure 3G in manuscript, changed for SFN

                                    )+
                            scale_y_continuous(breaks=c(0,5,10,20,40))+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none",
                                   axis.text.x=element_text(colour="white", size=scalefactor*34, family="sans"), #angle=45, hjust=1),
                                   axis.title.x=element_text(colour="white", size=scalefactor*38, family="sans"),
                                   axis.title.y=element_text(colour="black", size=scalefactor*38, family="sans"),
                                   #plot.title=element_text(colour="white",size=scalefactor*38, family="sans"),
                                   plot.subtitle=element_text(colour="white",size=scalefactor*38, family="sans", hjust=0.5),
                                   axis.line.x = element_blank(),
                                   axis.ticks.x = element_blank(), 
                 
                                   )

                     save_plot(filename = "unif_model_N", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                     N_sites<<- releaseProbPlot


  





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

                    #print(check_fitted)


                    print(N_sites_reg_model)

                    good_N_sites <- good_N_sites %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key)) %>%  dplyr::filter(trackROI_key %in% unique(fittable_ROIs$trackROI_key))#%>% dplyr::filter(N_min >= 0 )
                    generate_fit_predict <- generate_fit_predict %>% dplyr::filter(trackROI_key %in% unique(keep_ROIs$trackROI_key)) %>%  dplyr::filter(trackROI_key %in% unique(fittable_ROIs$trackROI_key))





                     ROI_list<- c("GluSnFR3-dish05-plate08-region02-ROI0018","GluSnFR3-dish05-plate09-region02-ROI0001")
                     #ROI_list<- unique(good_N_sites$trackROI_key) 
                     for (i in 1:length(ROI_list)) {


                            current_ROI = ROI_list[i]#"GluSnFR3-dish05-plate09-region02-ROI0002"
                                    print(paste0("loop iter is :", current_ROI) )

                            tmp_MV_score <- MV_per_synapse_df %>% dplyr::filter(trackROI_key == current_ROI)

                            tmp_sum_sq_resid <- get_sum_sq_resid %>% dplyr::filter(trackROI_key == current_ROI)
                            tmp_fit_predict<- generate_fit_predict %>% dplyr::filter(trackROI_key == current_ROI)
                            tmp_N_sites <- good_N_sites %>% dplyr::filter(trackROI_key == current_ROI) 
                            print(tmp_N_sites)
                            tmp_fit_values<- check_good_N_sites %>% dplyr::filter(trackROI_key == current_ROI)
                            unif_model = TRUE
                            nonunif_model = FALSE
                            model_score = "Uniform Model"
                            if(add_second_model == TRUE){
                            tmp_nonunif_fit_predict<- generate_nonunif_fit_predict %>% dplyr::filter(trackROI_key == current_ROI)
                            tmp_nonunif_fit_values<- nonunif_N_sites %>% dplyr::filter(trackROI_key == current_ROI)
                            print(tmp_nonunif_fit_predict)
                            print(tmp_nonunif_fit_values)
                            tmp_sum_sq_resid_nonunif <- get_sum_sq_resid_nonunif  %>% dplyr::filter(trackROI_key == current_ROI)

                             if(length(unique(tmp_sum_sq_resid$sum_sq_resid)) > 0 & length(unique(tmp_sum_sq_resid_nonunif$sum_sq_resid)) > 0){
                                          if(unique(tmp_sum_sq_resid$sum_sq_resid) <= unique(tmp_sum_sq_resid_nonunif$sum_sq_resid)){
                                            model_score = "Best Fit = Unif"
                                            unif_model = TRUE
                                            nonunif_model = TRUE
                                          } else {
                                            model_score = "Best Fit = NonUnif"
                                            unif_model = TRUE
                                            nonunif_model = TRUE
                                          }
                                      } else {
                                        model_score = "IDK what happened"
                                        unif_model = TRUE
                                        nonunif_model = TRUE
                                      }

                          }

                           
                                    if(max(unique(tmp_N_sites$mean_quanta), na.rm=TRUE) < 10){
                                          x_limit = 5
                                          #y_limit = 5
                                          break_vec = seq(0, x_limit, by = 2.5)
                                        
                                    } else {
                                          x_limit = 30
                                          #y_limit = 5
                                          break_vec = seq(0, x_limit, by = 5)
                                    }

                                    if(max(unique(tmp_fit_predict$var_quanta), na.rm=TRUE) < 5){
                            
                                            y_limit = 1.5
                                            #y_limit = 2*max(unique(tmp_fit_predict$var_quanta),na.rm=TRUE)
                                            y_break_vec = seq(0, 1.5, by = 0.5)
                                    }else{
                                      y_limit = 10
                                      y_break_vec = seq(0,y_limit, by = 2.5)
                                    }

                                    #print("checking whether models exist")

                                   

                            ###get values for P_rw at each Ca for printing on the graph
                            P05 = tmp_N_sites %>% dplyr::filter(Ca == "0pt5Ca")
                            P1 = tmp_N_sites %>% dplyr::filter(Ca == "1Ca")
                            P2 = tmp_N_sites %>% dplyr::filter(Ca == "2Ca")
                            P4 = tmp_N_sites %>% dplyr::filter(Ca == "4Ca")
                            print(P1)

                            numeric_P1 = round(unique(P1$P_rw),2 )# paste0("P[r_1Ca] ", round(unique(P1$P_rw),2 ))
                            numeric_P4 = round(unique(P4$P_rw),2 )#paste0("P[r_4Ca] ", round(unique(P4$P_rw),2 ))
                            numeric_Q = round(unique(tmp_N_sites$Q_w),2 )#paste0("Q ", round(unique(tmp_N_sites$Q_w),2 ))
                            numeric_N = round(unique(tmp_N_sites$N),2 )#paste0("N[sites] ", round(unique(tmp_N_sites$N),2 ))

                            label_P1 = bquote(P[r_1Ca] == .(numeric_P1))
                            label_P4 = bquote(P[r_4Ca] == .(numeric_P4))
                            label_Q = bquote(Q == .(numeric_Q))
                            label_N = bquote(N[sites] ==.(numeric_N))

                            N_sites<-ggplot(tmp_fit_predict, aes(x=mean_quanta, y=var_quanta, group=trackROI_key))+ 
                                                                        
                                                        {if(unif_model)geom_line( colour="black", size=1, alpha=0.7, lty="solid")}+
                                                        {if(nonunif_model == TRUE & add_second_model == TRUE)geom_line(data= tmp_nonunif_fit_predict, aes(x= mean_quanta, y=var_quanta, group=trackROI_key), colour="blue", size=0.8, alpha=0.7, lty="solid")}+
                                                        
                                                       
                                                        geom_point(data=tmp_N_sites,aes(x=mean_quanta,y=var_quanta, colour=Ca), size=8, alpha=0.9)+
                                                        
                                                        
                                                        #unif
                                                        {if(unif_model)geom_text(label = deparse(label_P1),  parse=TRUE, size=8, x= x_limit*0.7, y=y_limit*0.95)} + #parse=TRUE,
                                                        {if(unif_model)geom_text(label = deparse(label_P4), parse=TRUE, size=8, x= x_limit*0.7, y=y_limit*0.85)} +
                                                        {if(unif_model)geom_text(label = deparse(label_Q), parse=TRUE, size=8, x= x_limit*0.7, y=y_limit*0.75)} +
                                                        {if(unif_model)geom_text(label = deparse(label_N),  parse=TRUE, size=8, x= x_limit*0.7, y=y_limit*0.65)} +
                                                        #{if(unif_model)geom_text(label = paste0("resid.^2 = ", round(unique(tmp_sum_sq_resid$sum_sq_resid),4 )), size=7, x= x_limit*0.7, y=y_limit*0.55,colour='red')} +
                                                        #{if(nonunif_model == TRUE & add_second_model == TRUE)geom_text(label = paste0("Q = ", round(unique(tmp_nonunif_fit_values$Q),2 )), size=7, x= x_limit*0.7, y=y_limit*0.40, colour='blue')} +
                                                        #{if(nonunif_model == TRUE & add_second_model == TRUE)geom_text(label = paste0("N_sites = ", round(unique(tmp_nonunif_fit_values$N),2 )), size=7, x= x_limit*0.7, y=y_limit*0.30,colour='blue')} +
                                                        #{if(nonunif_model == TRUE & add_second_model == TRUE)geom_text(label = paste0("resid.^2 = ", round(unique(tmp_sum_sq_resid_nonunif$sum_sq_resid),4 )), size=7, x= x_limit*0.7, y=y_limit*0.20,colour='purple')} +
                                                        
                                                        
                                                        labs(x=expression(bar(x)~"("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"),
                                                                y=expression(sigma^2~"("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"),
                                                                tag="D",#"F", #Figure 3F in manuscript, changed for SFN poster 
                                                                subtitle=str_wrap("Single Bouton Mean-Variance Analysis",25) #model_score,
                                                                #subtitle=unique(tmp_MV_score$MV_score)
                                                                )+
                                                        coord_cartesian(xlim=c(0,x_limit),ylim=c(0,y_limit))+
                                                        
                                                        scale_colour_manual(values = color_override)+
                                                        scale_x_continuous(expand=c(0,0),breaks=break_vec)+
                                                        scale_y_continuous(expand=c(0,0))+     
                                                        
                                                        
                                                        theme_tufte()+
                                                        my.theme+
                                                        theme(legend.position = "none",
                                                                plot.margin = unit(c(1,1,1,1),"cm"),
                                                                )



                              save_plot(filename = paste0("mean-variance_nlsFits_", current_ROI), plot=N_sites, plot_width=plot_height*2,plot_height=plot_height*1.2, scale_f = scalefactor,dpi_val = 300)

                              if( current_ROI == "GluSnFR3-dish05-plate08-region02-ROI0018"){
                                         nam <- paste("N_sites_panelF", sep = "")
                                         assign(nam, N_sites,envir = .GlobalEnv)                                
                              }

                              }


                    wrap_val = 50
                   
                    #refactor
                    good_N_sites$groupfact <- factor(good_N_sites$Ca_mM)    
                    good_N_sites$Ca<- as.factor(as.character( (good_N_sites$Ca) ) )
                    levels(good_N_sites$Ca) <- new_levels

                    #disp_N_sites<- good_N_sites %>% dplyr::filter()    

                    releaseProbPlot<-ggplot(good_N_sites, aes(x=Ca_mM, y=P_rw))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=1.2,alpha=0.35)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            
                            stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white", alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression(P[r]),
                                    #title=expression("Binom. Uniform "*P[r]),
                                    tag="F"
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
                            coord_cartesian(ylim=c(0,1.0),xlim=c(0.25,4.25))+
                            
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")

                     save_plot(filename = "Pr_sina_plot", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                     Pr_distribution<<- releaseProbPlot

                     

                    #### RRP ratio??? #####
                     disp_N_sites$groupfact <- factor(disp_N_sites$Ca_mM)    
                      disp_N_sites$Ca<- as.factor(as.character( (disp_N_sites$Ca) ) )
                      levels(disp_N_sites$Ca) <- new_levels

                      disp_N_sites <- disp_N_sites %>% group_by(trackROI_key,Ca,Ca_mM) %>% mutate(RRP_ratio = mean_quanta/N)

                      bin=(4 - 0.5) / 50
                
                      releaseProbPlot<-ggplot(disp_N_sites, aes(x=Ca_mM, y=RRP_ratio))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=2,alpha=0.5)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            
                            stat_summary(geom="errorbar", fun.data=mean_se, width=bin*3, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.3, size=4,fill="white", alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression("("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"/N[sites]),
                                    subtitle=str_wrap("Estimated Proportion of RRP released per stimulus",25),
                                    tag="F" #"H" Figure 3H in manuscript, changed for SFN
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(0,0.5,1.0))+
                            coord_cartesian(ylim=c(0,1.2),xlim=c(0.25,4.25))+
                            
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")

                     save_plot(filename = "RRP_mobilized_plot", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                     RRP_proportion_distribution<<- releaseProbPlot


 N_sites_summary<- disp_N_sites %>% group_by(Ca, Ca_mM) %>% 
                                     summarise(n = n(),
                                               pop_mean_quanta = mean(mean_quanta),
                                               sd_quanta = sd(mean_quanta),
                                               se_quanta = sd_quanta/sqrt(n),

                                               mean_RRP_ratio = mean(RRP_ratio),
                                               sd_RRP_ratio = sd(RRP_ratio),
                                               se_RRP_ratio = sd_RRP_ratio/sqrt(n),
                                               
                                               mean_Pr = mean(P_rw),
                                               sd_Pr = sd(P_rw),
                                               se_Pr = sd_Pr/sqrt(n),

                                               mean_Qw = mean(Q_w),
                                               sd_Qw = sd(Q_w),
                                               se_Qw = sd_Qw/sqrt(n),

                                               mean_N_sites = mean(N),
                                               sd_N_sites = sd(N),
                                               se_N_sites = sd_N_sites/sqrt(n)
                                               )


#### potentially deprecated

                    #    light_good_N_sites <- good_N_sites %>% dplyr::filter(Ca_mM == 4)                        
                                                                                                
                    #    my_medians_N_sites <- light_good_N_sites %>%
                    #                                         group_by(Ca, Ca_mM) %>%
                    #                                       summarize(median_Q = median(Q_w,na.rm=TRUE),
                    #                                                 median_N_min = median(N_min,na.rm=TRUE) ) #,

                    # x_lower_bound = min(light_good_N_sites$N_min, na.rm=TRUE)
                    # x_upper_bound = max(light_good_N_sites$N_min, na.rm=TRUE)
                    # bin = (x_upper_bound-x_lower_bound)/25
                    # hist_plot <- ggplot(light_good_N_sites, aes(x = N_min)) +
                    #     geom_vline(data = my_medians_N_sites, aes(xintercept = median_N_min), size = 2, lty = "solid", colour="darkorchid2") +
                    #     geom_freqpoly(aes(y = ..count..), colour = "darkorchid2", binwidth = bin, alpha = 0.7, boundary = 0, size = 1.5) +
                    #     geom_histogram(aes(y = ..count..), fill = "darkorchid2", binwidth = bin, alpha = 0.4, boundary = 0, position = "identity", colour = "black", size = 0.8) +
                    #     #{if (quanta_vline) geom_vline(xintercept = 1, colour="red",lty="dashed",size=1,alpha=0.9 )}+
                    #     labs(
                    #       x = expression("Release sites, N"),
                    #       y = expression("Count"),
                    #       fill = "",
                    #       tag = "G",
                    #       #title="Binom. Release Sites, N"
                    #     ) +
                    #     # {if (color_switch == FALSE) scale_colour_brewer(palette = "Dark2")} +
                    #     # {if (color_switch) scale_colour_manual(labels = new_levels, values = color_override)} +
                    #     # {if (color_switch == FALSE) scale_fill_brewer(palette = "Dark2")} +
                    #     # {if (color_switch) scale_fill_manual(labels = new_levels, values = color_override)} +
                    #     coord_cartesian(ylim = c(0, 30*1.1), xlim = c(x_lower_bound, x_upper_bound)) +
                    #     scale_y_continuous(breaks = c(0, 10, 20,30), expand = c(0, 0)) +
                    #     theme_tufte() +
                    #     #facet_grid(Ca_expr ~ ., labeller = label_parsed) +
                    #     my.theme +
                    #     theme(
                    #       legend.position = "none",
                    #       legend.title = element_blank(),
                    #       strip.text.y.right = element_blank()
                    #     )

                    #   save_plot(
                    #     filename = paste0("histogram_N_sites", ""),
                    #     plot = hist_plot,
                    #     plot_width=plot_height*1.2,plot_height=plot_height*1.5, scale_f = scalefactor,dpi_val = 900)

                    #    N_sites_histogram<<- hist_plot


                    # ####    GET BINOMIAL Q #### ### compare to evoked iGlu 0.5 mM Ca2+

                    #      evoked_05Ca<- releaseProb_calc %>% dplyr::filter(Ca_mM == 0.5 ) %>% ungroup() %>% select(trackROI_key, mean_quanta) %>% mutate(expt = "Evoked")
                    #      names(evoked_05Ca)[names(evoked_05Ca) == 'mean_quanta'] <- 'Q_w'
                    #      Q_binom <- light_good_N_sites %>% ungroup() %>% select(trackROI_key, Q_w) %>% mutate(expt = "Binom. Fit")

                    #      Q_compare<- bind_rows(evoked_05Ca, Q_binom) %>% mutate(binary_expt = case_when(expt == "Evoked" ~ 0,
                    #                                                                                     expt == "Binom. Fit" ~ 1))

                    #      releaseProbPlot<-ggplot(Q_compare, aes(x=binary_expt, y=Q_w))+
                    #                             geom_sina(aes(group=binary_expt, colour=expt),size=1.2,alpha=0.35)+
                    #                             #stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                                                
                    #                             #geom_smooth(formula=y~x, method='loess',alpha=1,se=FALSE,size=3,colour="black")+              
                    #                             stat_summary(geom="errorbar", fun.data=mean_se, width=0.2, size=1, alpha=1, colour="black")+
                    #                             stat_summary(geom="point", fun.y=mean, shape=21,stroke=1.1, size=2.5,fill="white", alpha=1, colour="black")+
                                                
                    #                             labs( #x=expression("["*Ca^'2+'*"] (mM)"),
                    #                                     y=expression("Q ("*E[Delta*"F/F"]/S[Delta*"F/F"]*")"),
                    #                                     #title=expression(),
                    #                                     tag="H"
                    #                                     )+
                    #                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                    #                             {if(color_switch)scale_colour_manual(values=c("darkblue",color_override[1]))}+
                    #                             #guides(colour="none")+
                    #                             scale_x_continuous(breaks=c(0,1),labels=c("Evoked","Binom."))+
                    #                             #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
                    #                             #coord_cartesian(ylim=c(0,1.0),xlim=c(0.25,4.25))+
                    #                             #guides(colour = guide_legend(override.aes = list(colour = color_override,
                    #                              #                                                   size = c(4,4,4,4),
                    #                               #                                                  alpha=c(1,1,1,1)
                    #                                #                                                 )
                    #                                 #                                            )
                    #                                  #                                       )+                    
                    #                             theme_tufte()+
                    #                             my.theme+
                    #                             theme(legend.position = "none",
                    #                                   axis.title.x = element_text(colour="white", size=scalefactor*38, family="sans"))#,
                    #                                     #legend.text=element_text(colour="black", size=scalefactor*18, family="sans"))#,
                    #                             #               legend.justification = c('right','top'),
                    #                             #               axis.text.x=element_text(colour="black", size=20, family="sans"),
                    #                             #               axis.text.y=element_text(colour="black", size=20, family="sans"),
                    #                             #               axis.title=element_text(colour="black", size=24, family="sans"),
                    #                             #               axis.ticks.length=unit(.25, "cm"),
                    #                             #               axis.ticks = element_line(size=1),
                    #                             #               strip.text = element_text(colour="black", size = 16, family = "sans")
                    #                             #               )


                    #             #ggsave(filename=paste0("releaseProbPlot_v1"),plot=releaseProbPlot, device="jpeg",dpi=600, units="in",width=plot_height*1,height=plot_height*1)
                    #             save_plot(filename = "Q_sina_plot", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                    #                 Q_comparison<<- releaseProbPlot



MV_output<- list(MV_per_synapse_df = MV_per_synapse_df, good_N_sites = good_N_sites, fit_predict = generate_fit_predict, clean_N_sites=disp_N_sites, N_sites_summary = N_sites_summary)


}




#         current_ROI = good_ROIs[i]
#         tmp_sum_sq_resid <- get_sum_sq_resid %>% dplyr::filter(trackROI_key == current_ROI)
#         tmp_fit_predict<- generate_fit_predict %>% dplyr::filter(trackROI_key == current_ROI)
#         tmp_N_sites <- good_N_sites %>% dplyr::filter(trackROI_key == current_ROI) 
# #        print(tmp_N_sites)

#         tmp_fit_values<- check_good_N_sites %>% dplyr::filter(trackROI_key == current_ROI)
      
#        # if(unique(tmp_nonunif_fit_values$N) > 1000){
#        #     add_second_model = FALSE
#        #  } else {
#             add_second_model = FALSE
#        # }

#         if(max(unique(tmp_N_sites$mean_quanta), na.rm=TRUE) < 15){
#         x_limit = 15
            
#         } else {
#             x_limit = 25
#         }
#         y_limit = 2*max(unique(tmp_fit_predict$var_quanta),na.rm=TRUE)

#         #P1_nonunif = tmp_nonunif_fit_values %>% dplyr::filter(Ca == "1Ca")


