### plot_peakStats_activity_vs_marker.R
## written by Samuel T Barlow 12.6.23

# Goal of this script is to accept a summary data.frame, peakStats, and plot data visualizations for Figure 2 of the optical physiology manuscript. 
library(ggh4x)
library(ggpubr)
library(ggforce)
#library(Hmisc)
subDir <- "spont_and_evoked_stats_v3_just-in-case"   


source(paste0(path,"sub-functions/setFigureDirectory.R") )
source(paste0(path,"sub-functions/myTheme.R") )
#source(paste0(path,"sub-functions/def_stim_vlines.R") )
#source(paste0(path,"peakAnalysis/png_plotter_Figure2_v6.R") )
source(paste0(path,"peakAnalysis/read_synQuant_csv_v2.R") )
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R") )




plot_peakStats<- function(AP_stat_df = peak_stats, AP_df = df, spont_stat_df = spont_stats, spont_df = spont_df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops){



  						#### establish variables that will be used for generating datasets





                        scalefactor = 0.75
                        plot_height = 16
                        plot_width = 16#3/5*plot_height
                        lineplot_width = 16
                        lineplot_height = 24
                        scatterplot_width = 16
                        scatterplot_height = 16
                        color_switch = !is.null(color_override)

                        interFrame.var = as.vector(AP_df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)




###### Establish the dfs needed to plot averaged waveforms of evoked vs spontaneous at each Ca2+

 #### clean up dataframe by merging with peakStats

 						#initialize clean dataframes
 						 						remove_neurons = c("dish05-plate03-region02", "dish05-plate06-region01", "dish05-plate06-region02")#,"dish05-plate01-region01")
                        
                        AP_stat_df <- AP_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID),                                         
                        																		Ca_mM = case_when(Ca == "0pt5Ca" ~ 0.5,
                                                                                  Ca == "1Ca" ~ 1,
                                                                                  Ca == "2Ca" ~ 2,
                                                                                  Ca == "4Ca" ~ 4),
                                                                imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>%
                                                  dplyr::filter(!imaging_region %in% remove_neurons) #%>%
                                                  #dplyr::filter(tau_decay_ms < 100)
                        
                        
                        #### clean up dataframe by merging with peakStats

                        
                        light_df <- AP_df[,!(names(AP_df) %in% remove_col)]
                        
                       
                        ROIs_to_remove <- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are finite
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%
                                                 dplyr::filter( is.na(dFF))
                        ROIs_to_remove_timeskip<- light_df %>% dplyr::filter(stimEpoch == 1,absoluteTime < (firstStim + 0.5)) %>%        #only retain ROIs that are do not have a frame skip in time of interest
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%            
                                                 summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
                                                 dplyr::filter(maxInterframe > 0.18) %>%
                                                 select(ROI_key,maxInterframe) 



                        
                       ROIs_to_remove <- unique(ROIs_to_remove$ROI_key)     
                        ROIs_to_remove_timeskip<-unique(ROIs_to_remove_timeskip$ROI_key)                                             

                        clean_df <- light_df %>%  dplyr::filter(!ROI_key %in% ROIs_to_remove, !ROI_key %in% ROIs_to_remove_timeskip) %>% 
                                                  mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>%
                                                  dplyr::filter(!imaging_region %in% remove_neurons)


#### checking things		stat_ROI
												#check_shorter_taus<- AP_stat_df %>% dplyr::filter(tau_decay_ms <= 50)
                        stat_ROIs <- unique(AP_stat_df$ROI_key)
                        print("These should be the unique tau_decay_ms")
                        print(unique(AP_stat_df$tau_decay_ms))

##### check plot avgs early
                        plot_avgs_test <- AP_stat_df %>%
																			  	group_by(Ca, segmentation,protocol) %>%
																			  summarise(n = n(),
																			         mean_amplitude = mean(amplitude),
																			         sd_amplitude = sd(amplitude),
																			         se_amplitude = sd_amplitude / sqrt(n),
																			         
																			         mean_tau = mean(tau_decay_ms),
																			         sd_tau = sd(tau_decay_ms),
																			         se_tau = sd_tau / sqrt(n),

																			         mean_t_rise = mean(t_rise),
																			         sd_t_rise = sd(t_rise),
																			         se_t_rise = sd_t_rise / sqrt(n),

																			         mean_t_decay = mean(t_decay),
																			         sd_t_decay = sd(t_decay),
																			         se_t_decay = sd_t_decay / sqrt(n),

																			         mean_t_half = mean(t_half),
																			         sd_t_half = sd(t_half),
																			         se_t_half = sd_t_half / sqrt(n),

																			         mean_dt = mean(interSpike_ms,na.rm=TRUE),
																			         sd_dt = sd(interSpike_ms,na.rm=TRUE),
																			         se_dt = sd_dt / sqrt(n),
																			         )#,
									         
									         print(plot_avgs_test)


									           plot_avgs_test_fix_infinite <- AP_stat_df %>%
																			  	group_by(Ca, segmentation,protocol) %>%
																			  summarise(n = n(),
																			         mean_amplitude = mean(amplitude),
																			         sd_amplitude = sd(amplitude),
																			         se_amplitude = sd_amplitude / sqrt(n),
																			         
																			         mean_tau = mean(tau_decay_ms),
																			         sd_tau = sd(tau_decay_ms),
																			         se_tau = sd_tau / sqrt(n),

																			         mean_t_rise = mean(t_rise),
																			         sd_t_rise = sd(t_rise),
																			         se_t_rise = sd_t_rise / sqrt(n),

																			         mean_t_decay = mean(t_decay),
																			         sd_t_decay = sd(t_decay),
																			         se_t_decay = sd_t_decay / sqrt(n),

																			         mean_t_half = mean(t_half),
																			         sd_t_half = sd(t_half),
																			         se_t_half = sd_t_half / sqrt(n),

																			         mean_dt = mean(interSpike_ms[!is.infinite(interSpike_ms)], na.rm = TRUE),
																			         sd_dt = sd(interSpike_ms[!is.infinite(interSpike_ms)], na.rm = TRUE),
																			         se_dt = sd_dt / sqrt(n),
																			         )#,
									         




                        clean_peaks<- suppressMessages(clean_df %>% dplyr::filter(windowedPeakID != "NotPeak", ROI_key %in% stat_ROIs) )    #left_join(clean_df,AP_stat_df) %>% dplyr::filter(windowedPeakID != "NotPeak") ) 
                        timezero_df<- suppressMessages(clean_peaks %>%  group_by_at(groupers) %>% summarise(time_zero = unique(firstStim[!is.na(firstStim)] ) ) )
                        clean_peaks_normTime<- suppressMessages(left_join(clean_peaks,timezero_df) %>% mutate(new_normTime = absoluteTime - time_zero) ) #%>% dplyr::filter(ROI_key %in% ROIs_with_peaks))
     					

     					#get levels
						Ca_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
						protocol_levels = c("evoked","spontaneous")
						segmentation_levels = c("marker","activity")

						#convert data for new_levels
     					clean_peaks_normTime$Ca<- as.factor(as.character( (clean_peaks_normTime$Ca) ) )
     					AP_stat_df$Ca<- as.factor(as.character( (AP_stat_df$Ca) ) )
     					#assign new levels
                        levels(AP_stat_df$Ca) <- Ca_levels
                        levels(clean_peaks_normTime$Ca) <- Ca_levels
						#levels(AP_stat_df$protocol) <- Ca_levels
                        #levels(clean_peaks_normTime$protocol) <- Ca_levels
                        

                        spont_stat_df<- spont_stat_df %>% ungroup() %>% group_by(Ca) %>%
                        																								dplyr::filter(tau_decay_ms >= 10,
                        																															tau_decay_ms <= 200,
                        																															tau_decay_ms <= quantile(tau_decay_ms, 0.985, na.rm=TRUE), 
                                                                                			tau_decay_ms >= quantile(tau_decay_ms, 0.015, na.rm=TRUE),
                        																								 							amplitude >= 0.05,
                                                                                			t_half > 0,
                                                                                			t_rise > 0,
                                                                                			t_decay > 0)
												spont_ROIs_to_keep <- unique(spont_stat_df$ROI_key)




                        clean_spont_df <- data.frame(spont_df) %>% ungroup() %>% dplyr::filter(windowed_normTime >= -0.055) %>%
                        							   				 dplyr::filter(windowed_normTime <= 0.2) %>%
                        							   				 dplyr::filter(ROI_key %in% spont_ROIs_to_keep)

                        names(clean_peaks_normTime)[names(clean_peaks_normTime) == 'new_normTime'] <- 'windowed_normTime'

                    
                        


#### generate ROI count from raw data

						getROIs<- clean_df %>% dplyr::filter(!ROI_key %in% ROIs_to_remove) %>% group_by(vid_key, segmentation,Ca, Ca_mM) %>% summarise(ROI_count = length(unique(ROI_key)))
						getROIs$segmentation <- factor(getROIs$segmentation, levels = segmentation_levels)

						ROIs_with_peaks<- clean_df %>% dplyr::filter(!ROI_key %in% ROIs_to_remove,stimEpoch == 1) %>% group_by(segmentation,Ca, Ca_mM, ROI_key) %>% summarise(peaks_detected = length(unique(peakID))-1,
																																					peak_positive = case_when(peaks_detected > 0 ~ 1,
																																												peaks_detected == 0 ~ 0)
																																				)


						area_imaging_region_um = 25.6 * 25.6
						check_ROIs_per_vid <- suppressMessages(getROIs %>% group_by(segmentation, Ca) %>% mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) %>% 
																										  summarise(total_imaging_regions = length(unique(imaging_region)),
									                                                                                    total_vids = length(unique(vid_key)),
									                                                                                    sum_ROIs = sum(ROI_count),
									                                                                                    ROIs_per_vid = sum_ROIs/total_vids,
									                                                                                    mean_ROI_count_per_vid = mean(ROI_count),
									                                                                                    sd_ROI_count = sd(ROI_count),
									                                                                                    se_ROI_count = sd_ROI_count/sqrt(sum_ROIs),
									                                                                                    as_area_mean_ROI = mean_ROI_count_per_vid/area_imaging_region_um,
									                                                                                    as_area_sd_ROI = sd_ROI_count/area_imaging_region_um ) )

						print(getROIs)







##### Establish average waveform plotting
vline_df <- data.frame(protocol_fixed = c('Spontaneous','Evoked'), xintercept = c(NA,0))

try_normalize<- suppressMessages(clean_spont_df %>% group_by_at(groupers) %>% summarise(max_norm_time = max(windowed_normTime,na.rm=TRUE)) )
check<- suppressMessages(left_join(clean_spont_df, try_normalize) %>% 
							group_by_at(groupers) %>% 
							dplyr::filter(windowed_normTime >= (unique(max_norm_time) - 0.05), windowed_normTime <= unique(max_norm_time))  %>% 
							summarise(mean_offset = mean(dFF,na.rm=TRUE))
							)
clean_spont_df <- suppressMessages(left_join(clean_spont_df,check) %>%
									group_by_at(groupers) %>% 
									dplyr::filter(!is.na(peakgrp) ) %>%
									mutate(offset_dFF = dFF - mean_offset)
									 )


try_normalize<- suppressMessages(clean_peaks_normTime %>% group_by_at(groupers) %>% summarise(max_norm_time = max(windowed_normTime,na.rm=TRUE)) )
check<- suppressMessages(left_join(clean_peaks_normTime, try_normalize) %>% 
							group_by_at(groupers) %>% 
							dplyr::filter(windowed_normTime >= (unique(max_norm_time) - 0.05), windowed_normTime <= unique(max_norm_time))  %>% 
							summarise(mean_offset = mean(dFF,na.rm=TRUE))
							)
clean_peaks_normTime <- suppressMessages(left_join(clean_peaks_normTime,check) %>%
									group_by_at(groupers) %>%
									dplyr::filter(windowed_normTime >= -0.055, windowed_normTime <= 0.4) %>%  
									#dplyr::filter(!is.na(amplitude) ) %>%
									mutate(offset_dFF = dFF - mean_offset)
									 )


hline_df <- data.frame(protocol_fixed = c('Spontaneous','Evoked','Spontaneous','Evoked','Spontaneous','Evoked','Spontaneous','Evoked'),
												Ca = c("\"0.5 mM \" * Ca^\"2+\"","\"0.5 mM \" * Ca^\"2+\"",
																 "\"1 mM \" * Ca^\"2+\"", "\"1 mM \" * Ca^\"2+\"",
																 "\"2 mM \" * Ca^\"2+\"", "\"2 mM \" * Ca^\"2+\"",
																 "\"4 mM \" * Ca^\"2+\"","\"4 mM \" * Ca^\"2+\""),
												vline_val = c(0.03555, 0.05536,
																			0.03657, 0.07678,
																			0.03869, 0.10081,
																			0.03337, 0.12883 ###avg t_decay
																			)
												)


 span_max = 0.55 - (-0.055)
 span_factor = span_max/interFrame.var
 span_arg = 1/span_factor
merge_df = bind_rows(clean_peaks_normTime,clean_spont_df) %>% mutate(protocol_fixed = case_when(protocol == "singleAP" ~ "Evoked",
																																																protocol == "spont" ~ "Spontaneous")) %>% 
																															dplyr::filter(windowedPeakID != "NotPeak") %>%
																															dplyr::filter()
																															

																																				#offset_dFF	
						avgtrace_Plot<-ggplot(merge_df, aes(x=windowed_normTime, y = dFF,  group=interaction(Ca,protocol_fixed), colour=Ca))+
						                            #geom_path(colour="grey61",size=0.5,alpha=0.4)+
						                            #geom_vline(xintercept = c(0.015,0.028,0.04,0.06),colour='black',lty='dashed')+
						                            #geom_vline(data=hline_df, aes(xintercept = vline_val, group=protocol_fixed), lty="solid",colour='red',size=1)+
						                            #geom_smooth(aes(group=Ca), span=span_arg)+
						                            #stat_summary_bin(aes(group=Ca), geom = "errorbar", fun.data=mean_se, binwidth=interFrame.var,size=0.2,alpha=1,colour="grey11")+
						                            stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.5,alpha=1)+
						                            
						                            #stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.1,alpha=0.6)+
						                            geom_vline(data=vline_df, aes(xintercept = xintercept, group=protocol_fixed), lty="dashed",colour='grey21',size=0.9)+
						                            
						                           	labs( x="Normalized time (s)",
						                                  y=expression(Delta*"F/F"),
						                                  subtitle="Averaged iGluSnFR3 Responses",
						                                  #colour="Segmented by:",
					                                                    tag = "A")+
						                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                            {if(color_switch)scale_colour_manual(labels=Ca_levels, values=color_override)}+
                                                                #scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						                            scale_x_continuous(labels = c("0", "0.5"),breaks = c(0,0.5),limits=c(-0.055,0.65))+
						                            theme_tufte()+
						                            #guides(colour="none",fill="none")+
						                            coord_cartesian(ylim=c(-0.01,3.5),xlim=c(-0.055,0.55))+
						                            my.theme+
						                            #facet_grid(~Ca + protocol, labeller = label_parsed)+
						                            facet_nested_wrap(~Ca + protocol_fixed, nrow = 1, ncol=8, labeller = labeller(Ca = function(x) {rep("", length(x))}))+
						                            theme(legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
						                            	legend.position=c(0.2,0.7),
						                            	legend.justification= "right",
						                            	panel.spacing = unit(0.25, "cm"),
						                            	strip.clip = "off")#strip.text.x = element_text(colour="black", size=scalefactor*28, family="sans"))

						                           avg_peaks  <<- avgtrace_Plot
						                            #       rm(avg_PP_tracePlot
						                                

						                             save_plot(filename=paste0("averagePeaks_"),plot=avgtrace_Plot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

### check omit long taus



                       
### stats as boxplots




Ca_str = unique(AP_stat_df$Ca)

merge_stats<- bind_rows(AP_stat_df,spont_stat_df) %>%  mutate(protocol_fixed = case_when(protocol == "singleAP" ~ "Evoked",
																								protocol == "spont" ~ "Spontaneous"),
																Ca_mM = case_when(Ca == "\"0.5 mM \" * Ca^\"2+\"" ~ 0.5,
																					Ca == "\"1 mM \" * Ca^\"2+\"" ~  1,
																					Ca == "\"2 mM \" * Ca^\"2+\"" ~ 2,
																					Ca == "\"4 mM \" * Ca^\"2+\"" ~ 4))



                        plot_avgs <- merge_stats %>%
									  	group_by(Ca, segmentation, protocol, protocol_fixed) %>%
									  summarise(n = n(),
									         mean_amplitude = mean(amplitude),
									         sd_amplitude = sd(amplitude),
									         se_amplitude = sd_amplitude / sqrt(n),
									         
									         mean_tau = mean(tau_decay_ms),
									         sd_tau = sd(tau_decay_ms),
									         se_tau = sd_tau / sqrt(n),

									         mean_t_rise = mean(t_rise),
									         sd_t_rise = sd(t_rise),
									         se_t_rise = sd_t_rise / sqrt(n),

									         mean_t_decay = mean(t_decay),
									         sd_t_decay = sd(t_decay),
									         se_t_decay = sd_t_decay / sqrt(n),

									         mean_t_half = mean(t_half),
									         sd_t_half = sd(t_half),
									         se_t_half = sd_t_half / sqrt(n),

									         mean_dt = mean(interSpike_ms),
									         sd_dt = sd(interSpike_ms),
									         se_dt = sd_dt / sqrt(n),
									         )#,
									         

print(unique(merge_stats$Ca))


####


#### PICKING UP HERE, 4.29.24 STB ####


###
 light_spont_avg <- merge_stats %>% group_by(Ca, Ca_mM) %>%
 																		dplyr::filter(protocol_fixed == "Spontaneous") %>%
 																		summarise(spont_amplitude = mean(amplitude,na.rm=TRUE) ) %>%
                                                    ungroup() %>%
                                                    select(Ca_mM,spont_amplitude)
light_evoked_stat_df <- merge_stats %>% dplyr::filter(protocol_fixed == "Evoked")
quanta_stat_df <- suppressMessages(left_join(light_evoked_stat_df, light_spont_avg) %>% 
                                                                            mutate(quanta_calc = amplitude/spont_amplitude) )

  mean_quanta_stat_df <- suppressMessages( quanta_stat_df %>% 
                                                         ungroup() %>%
                                                         group_by(Ca,Ca_mM,vid_key,trackROI_key,ROINumber) %>%
                                                         summarise(mean_quanta = mean(quanta_calc, na.rm=TRUE) )
                                                         ) 
  							mean_quanta_stat_df$groupfact <- factor(mean_quanta_stat_df$Ca_mM)

								bin=(4 - 0.5)	/ 50
                releaseProbPlot<-ggplot(mean_quanta_stat_df, aes(x=Ca_mM, y=mean_quanta))+
                            geom_sina(aes(group=groupfact, colour=groupfact),size=2,alpha=0.5)+
                            stat_summary(geom="line", fun.y=mean, size=1, alpha=0.85, colour="black")+
                            
                            geom_hline(yintercept = 1, size=0.7,lty="dashed",colour="black")+
                            #geom_smooth(formula=y~x, method='loess',alpha=1,se=FALSE,size=3,colour="black")+              
                            stat_summary(geom="errorbar", fun.data=mean_se, width=bin*3, size=1, alpha=1, colour="black")+
                            stat_summary(geom="point", fun.y=mean, shape=21,size=4,stroke=1.3,fill="white",,alpha=1, colour="black")+
                            
                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                    y=expression("E"[Delta*"F/F"] / "S"[Delta*"F/F"]),
                                    subtitle = str_wrap("Estimated SVs released per stimulus",25),
                                    tag="C" #"E" Figure 3E in manuscript, changed for SFN
                                    )+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=Ca_levels, values=color_override)}+
                            #guides(colour="none")+
                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            scale_y_continuous(breaks=c(1,10,20,30))+
                            coord_cartesian(ylim=c(0,35),xlim=c(0.25,4.25))+
                                                
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none",
                            			axis.text.x = element_text(colour="black", size=scalefactor*34, family="sans"))#, angle=30,hjust=1))#,
                            #               legend.justification = c('right','top'),
                            #               axis.text.x=element_text(colour="black", size=20, family="sans"),
                            #               axis.text.y=element_text(colour="black", size=20, family="sans"),
                            #               axis.title=element_text(colour="black", size=24, family="sans"),
                            #               axis.ticks.length=unit(.25, "cm"),
                            #               axis.ticks = element_line(size=1),
                            #               strip.text = element_text(colour="black", size = 16, family = "sans")
                            #               )

            save_plot(filename = "quanta_plot_v1", plot=releaseProbPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
            
            #ggsave(filename=paste0("quanta_plot_v1"),plot=releaseProbPlot, device="jpeg",dpi=600, units="in",width=plot_height*1,height=plot_height*1)
            
            quanta_conv<<- releaseProbPlot
                                                

# Set global parameters
#lower_bound_x <- 0
#upper_bound_x <- 10
#bin <- (upper_bound_x - lower_bound_x) / 50
new_color_override <- c("black", "darkorange3")

               my_medians <- merge_stats %>%
									  	group_by(Ca, protocol_fixed) %>%
									  summarise(median_amplitude = median(amplitude),
									  			median_tau = median(tau_decay_ms),
									  			median_interSpike =median(interSpike_ms,na.rm=TRUE),
									  			median_t_rise = median(t_rise,na.rm=TRUE),
											  			median_t_half = median(t_half,na.rm=TRUE),
											  			median_t_decay = median(t_decay,na.rm=TRUE) ) #,



tag_var = c("B","C","D","E","F")
# HISTOGRAM FUNCTION
create_histogram <- function(data, median_data, x_var, lower_bound_x,upper_bound_x, median_var, title, filename,global_output,tmp_tag_var=NULL,legend_bool = TRUE ) {


	theme_histo = function(bool = TRUE) {
					       if(bool == TRUE) {
					       	theme(	legend.position = c(0.6,0.95),
					       			legend.title = element_blank(),
					       			strip.text.y.right = element_blank()
					       			)
					       	} else {
					       		theme(legend.position = "none",
					       				legend.title = element_blank(),
					       				strip.text.y.right = element_blank()
					       			) 

					       	} 

					       
					}

  print(paste0("creating histogram for :", x_var) )

  bin = (upper_bound_x - lower_bound_x)	/ 50
  histo <- ggplot(data, aes_string(x = x_var, group = "protocol_fixed", fill = "protocol_fixed", colour = "protocol_fixed")) +
    geom_vline(data = median_data, aes_string(xintercept = median_var, colour = "protocol_fixed"), size = scalefactor * 2, lty = "solid",show_guide = FALSE,alpha=0.9) +
    geom_freqpoly(aes(y = ..ncount..), binwidth = bin, alpha = 0.8, boundary = 0, size = scalefactor * 1,show_guide = FALSE) +
    geom_histogram(aes(y = ..ncount.., alpha = protocol_fixed), binwidth = bin, boundary = 0, position = "identity", colour = "black", size = scalefactor * 0.5,show_guide = TRUE) +
    #geom_vline(xintercept = 15,lty="dashed",colour="black")+ ### checking distribution cutoffs
    
    {if(legend_bool)geom_text(aes(label = Ca), x = upper_bound_x+(0.05*upper_bound_x), y = 0.45, size=8,hjust=1,check_overlap = TRUE,parse=TRUE,show_guide=FALSE)}+
    labs(x = title,
         y = expression(N / N[max]),
		 tag = tmp_tag_var
         ) +
    {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
	{if(color_switch)scale_colour_manual( values=new_color_override)}+
    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
	{if(color_switch)scale_fill_manual( values=new_color_override)}+
	scale_alpha_manual(values = c("Evoked"= 0.4,"Spontaneous"=0.6))+
    coord_cartesian(ylim = c(0, 1.3), xlim = c(lower_bound_x, upper_bound_x)) +
    scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
    theme_tufte() +
    facet_grid(Ca ~ ., labeller = label_parsed) +
    my.theme +
    theme_histo(bool = legend_bool)

    # theme(legend.position = "none",
    #       legend.title = element_blank(),
    #       strip.text.y.right = element_blank())

  save_plot(filename = paste0(filename, ""), plot = histo, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)


  	if(global_output ==TRUE){
				nam <- filename
                assign(nam, histo,envir = .GlobalEnv)
            }

}

output_graphs=TRUE
# # Create histograms
create_histogram(merge_stats, my_medians, "amplitude", 0,10, "median_amplitude", expression("Peak " * Delta * "F/F"), "histo_amplitude",output_graphs,tag_var[2],TRUE)
create_histogram(merge_stats, my_medians,"tau_decay_ms",0,200, "median_tau", expression(tau[decay] * " (ms)"), "histo_tau_decay_ms",output_graphs,tag_var[3],FALSE)
create_histogram(merge_stats, my_medians,"t_rise",0,40, "median_t_rise", expression(t[rise] * " (ms)"), "histo_t_rise",FALSE,NULL,FALSE)
create_histogram(merge_stats, my_medians,"t_decay",0,400, "median_t_decay", expression(t[decay] * " (ms)"), "histo_t_decay",FALSE,NULL,FALSE)
create_histogram(merge_stats, my_medians,"t_half",0,150, "median_t_half", expression(t[1/2] * " (ms)"), "histo_t_half",output_graphs,tag_var[4],FALSE)
      


##### SCATTER PLOTS

create_scatterplot <- function(data, y_var, my_breaks,lower_bound_x,upper_bound_x, title, filename,global_output,tmp_tag_var = NULL,legend_bool=TRUE ) {



					if(legend_bool == TRUE) {
								subtitle = "Average Amplitude"
					} else {
						subtitle = NULL
					}

					theme_scatter = function(bool = TRUE) {
					       if(bool == TRUE) {
					       	theme(axis.text.x = element_text(colour="black", size=scalefactor*34, family="sans"),
					       			legend.position = c(0.4,0.9)
					       			)
					       			
					       	} else {
					       		theme(axis.text.x = element_text(colour="black", size=scalefactor*34, family="sans"),
					       			legend.position = "none"
					       			) 
					       			#subtitle = ""
					       	} 

					       
					}

					scalefactor=0.75
				  print(paste0("creating scatterplot for :", y_var) )
				  bin = (4 - 0.5)	/ 50
				  ymin = first(my_breaks)
				  ymax = last(my_breaks)  
				  label.p = 0.85*ymax
				  #legend.pos = ifelse(legend_bool == TRUE, c(0.4,0.95),"none")
				  #print(legend.pos)
				  data <- data %>% mutate(char_Ca = case_when(Ca_mM == 0.5 ~ "0.5",
				  												Ca_mM == 1 ~ "1",
				  												Ca_mM == 2 ~ "2",
				  												Ca_mM == 4 ~ "4"))
				  scatter<-ggplot(data, aes_string(x="Ca_mM", y=y_var,group="protocol_fixed",colour="protocol_fixed"))+
				                            stat_summary(geom="errorbar", fun.data=mean_se, width=bin*3, size=1, alpha=1)+
				                            stat_summary(geom="line", fun.y=mean, size=1, alpha=1)+
				                            stat_summary(geom="point", fun.y=mean, size=4, alpha=1, stroke=1.3,shape = 21, fill = "white")+
				                            #stat_compare_means(aes_string(group= "protocol_fixed"), label = "p.signif", method = "wilcox.test",label.y=label.p,size=8,show.legend=FALSE)+
				                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
																					y=title,
																					subtitle = subtitle,
				                                    tag = tmp_tag_var)+
				                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
				                            {if(color_switch)scale_colour_manual(values=new_color_override)}+
				                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+#, labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                                            scale_y_continuous(breaks=my_breaks)+
				                            coord_cartesian(ylim=c(ymin,ymax ),xlim=c(0.5,4))+
				                            theme_tufte()+
				                            my.theme+
				                            theme_scatter(bool = legend_bool)

				            save_plot(filename = paste0(filename, ""), plot = scatter, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)




				            if(global_output ==TRUE){
								nam <- filename
				                assign(nam, scatter,envir = .GlobalEnv)
				            }
}

break_vec = c(0,1,2,3,4)
create_scatterplot(merge_stats, "amplitude", break_vec,0,4, expression("Peak " * Delta * "F/F"), "nostats_scatter_amplitude",TRUE,"B",TRUE) #tag_var[1]
break_vec = c(25,50,75,100)
create_scatterplot(merge_stats,"tau_decay_ms",break_vec,0,200, expression(tau[decay] * " (ms)"), "nostats_scatter_tau_decay",TRUE,tag_var[4],FALSE)
break_vec = c(10,15,20,25,30,35)
create_scatterplot(merge_stats,"t_rise",break_vec,0,50,  expression(t[rise] * " (ms)"), "nostats_scatter_t_rise",TRUE,tag_var[3], FALSE)
break_vec = c(50,100,150,200,250)
create_scatterplot(merge_stats,"t_decay",break_vec,0,400, expression(t[decay] * " (ms)"), "nostats_scatter_t_decay",TRUE,tag_var[5],FALSE)
break_vec = c(20,40,60,80)
create_scatterplot(merge_stats,"t_half",break_vec,0,200,  expression(t[1/2] * " (ms)"), "nostats_scatter_t_half",TRUE,tag_var[5],FALSE)




# # ##### GET SYNQUANT CORRELATIONS

# AP_SQ_dir = "Y:\\Sam/paper1_datasets/singleAP_v2/imageFiles/zstack/zstackOutputs_sliced_SQ"
# spont_SQ_dir = "Y:\\Sam/paper1_datasets/spont_v1/imageFiles/zstack/zstackOutputs_SQ"

# file_pattern = "Ch2Results.csv"


# AP_SQ_data<- SQ_reader(AP_SQ_dir,file_pattern,"Evoked")
# spont_SQ_data<- SQ_reader(spont_SQ_dir,file_pattern,"Spontaneous")
# print(head(spont_SQ_data))

# #merge_AP_SQ_stats<- suppressMessages(left_join(AP_stat_df,AP_SQ_data))
# #merge_spont_SQ_stats<- suppressMessages(left_join(spont_stat_df,spont_SQ_data))
# #print(head(merge_spont_SQ_stats))

# check_spont_stats_freq<- merge_stats %>% dplyr::filter(protocol_fixed == "Spontaneous") %>% group_by(dish,region,exposeNum,vid_key,Ca,protocol_fixed,ROI_key,ROINumber) %>% 
# 										   summarise(peak_freq = length(unique(peakgrp))/31)# %>% 
# 										   #ungroup() %>% 
# 										   #group_by(dish,region,Ca,protocol_fixed,ROINumber) %>%
# 										   #summarise(overall_peak_freq = mean(peak_freq,na.rm=TRUE))
# merge_spont_SQ_freq<- left_join(check_spont_stats_freq, spont_SQ_data) %>% 
# 											dplyr::filter(!is.na(protocol_fixed)) %>% 
# 											group_by(dish,region,protocol_fixed,Ca) %>% 
# 											mutate(mean_IntDen = mean(IntDen, na.rm=TRUE),
# 													sd_IntDen = sd(IntDen,na.rm=TRUE),
# 													z_IntDen = (IntDen - mean_IntDen)/sd_IntDen 
# 													)


# # merge_SQ_stats<- bind_rows(merge_AP_SQ_stats,merge_spont_SQ_stats) %>% dplyr::filter(!is.na(protocol_fixed)) %>% group_by(dish,region,protocol_fixed,Ca) %>% mutate(mean_IntDen = mean(IntDen, na.rm=TRUE),
# # 																																			  					sd_IntDen = sd(IntDen,na.rm=TRUE),
# # 																																			  					z_IntDen = (IntDen - mean_IntDen)/sd_IntDen 
# # 																																			  	)
# rm(merge_AP_SQ_stats,merge_spont_SQ_stats)

# #print(names(merge_SQ_stats))
# source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 


# create_correlation_plot <- function(data, x_var, y_var, x_title,y_title,  filename) {




# 					scalefactor=0.75
# 				  print(paste0("creating correlation scatterplot for : x=", x_var, " and y=", y_var) )
# 				  #bin = (4 - 0.5)	/ 50
# 				  #ymin = first(my_breaks)
# 				  #ymax = last(my_breaks)  
# 				  #label.p = 0.85*ymax
# 				  #legend.pos = ifelse(legend_bool == TRUE, c(0.4,0.95),"none")
# 				  #print(legend.pos)
# 				  #data <- data %>% mutate(char_Ca = case_when(Ca_mM == 0.5 ~ "0.5",
# 				 # 												Ca_mM == 1 ~ "1",
# 				 # 												Ca_mM == 2 ~ "2",
# 				 # 												Ca_mM == 4 ~ "4"))
# 				 # 
# 				 corr_plot<-ggscatter(data, x = x_var, y = y_var, 
# 									          add = "reg.line", conf.int = TRUE, 
# 									          cor.coef = FALSE, cor.method = "pearson",
# 									          size = 1,
# 									          alpha=0.4,
# 									          #xlab = x_title, ylab = y_title, 
# 									          color="protocol_fixed",
# 									          palette = c(Evoked = "black", Spontaneous = "darkorange3")
# 									          )+

#                         		ggpubr::stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), color = "black", geom = "label") +
# 	                            labs( x=x_title,
# 	                                    y=y_title
# 	                                    )+
# 	                            scale_x_continuous(breaks = c(0,2,4,6,8))+
# 	                            #{if(color_switch==FALSE)scale_colour_brewer(palette="Dark2")}+
# 	                            #{if(color_switch)scale_colour_manual(values=new_color_override)}+
# 	                            theme_tufte()+
# 	                            my.theme+
# 	                            facet_grid(.~Ca, labeller=label_parsed)
# 	        					ggsave(filename = paste0(filename, ""), plot = corr_plot, device = "jpeg", dpi = 600, bg = "white", units = "in", width = plot_width*scalefactor, height = plot_height*scalefactor)
                   
# }

# #create_correlation_plot(merge_SQ_stats, "z_IntDen","amplitude",expression("Integrated Density z-score"), expression("Peak " * Delta * "F/F"),"corr_plot_amp_vs_IntDen_v3")
# create_correlation_plot(merge_spont_SQ_freq, "z_IntDen","peak_freq",expression("Integrated Density z-score"), expression("Mini frequency (Hz"^-1*")"),"corr_plot_minifreq_vs_IntDen_v1")
			










list.output = list(ROI_totals = getROIs, ROI_probability = ROIs_with_peaks, ROIs_per_vid_summary = check_ROIs_per_vid, Medians = my_medians, Averages = plot_avgs, clean_stats =merge_stats,clean_df_normTime = merge_df, test_avgs = plot_avgs_test_fix_infinite)
list.output
#vline_df
}









# ### as boxplots
#     				boxplot_amplitude<-ggplot(merge_stats, aes(x=Ca_mM, y=amplitude, group=interaction(Ca,protocol_fixed), colour=Ca))+#linetype=segmentation,
#                                              #stat_ecdf(size=2,alpha=1, lineend="round")+
#                                              #geom_sina(alpha=0.3)+
#                                              geom_boxplot(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75),size=2)+
#                                             labs( y=expression("Peak "*Delta*"F/F"),
#                                                     x=expression(Ca^"2+"*" (mM)"),#"Cumulative Density",
#                                                     #title="Peak amplitude",
#                                                     #colour="Segmented by:"
#                                                     )+
                                            
#                                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                             {if(color_switch)scale_colour_manual(labels=Ca_levels, values=color_override)}+
#                                             coord_cartesian(ylim=c(0,10))+#,xlim=c(0,5))+
#                                             #scale_y_continuous(breaks=c(0,0.5,1.0))+
#                                             #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
#                                             theme_tufte()+
#                                             facet_grid(~protocol_fixed, labeller = label_parsed)+
#                                             my.theme#+
#                                             #theme(legend.title = element_text(colour="black", size=28, family="sans"),
#                                             	  #strip.text.y.right = element_text(angle = 0),
#                                             #	  axis.text.x=element_blank())
#                                             #legend.position = c(.95, .95),
#                                                   #legend.justification = c("right", "bottom"))#,
#                                                    # legend.margin = margin(6, 6, 6, 6))
                       
#                         ggsave(filename=paste0("amplitude_boxplot"),plot=boxplot_amplitude, device="jpeg",dpi=600, bg="white",units="in",width=plot_height,height=plot_width)

#                     boxplot_tau_decay<-ggplot(merge_stats, aes(x=Ca_mM, y=tau_decay_ms, group=interaction(Ca,protocol_fixed), colour=Ca))+#linetype=segmentation,
#                                              #stat_ecdf(size=2,alpha=1, lineend="round")+
#                                              #geom_sina(alpha=0.3)+
#                                              geom_boxplot(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75),size=2)+
#                                             labs( y=expression(tau[decay]*" (ms)"),
#                                                     x="",#"Cumulative Density",
#                                                     #title="Peak amplitude",
#                                                     #colour="Segmented by:"
#                                                     )+
                                            
#                                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                                             {if(color_switch)scale_colour_manual(labels=Ca_levels, values=color_override)}+
#                                             coord_cartesian(ylim=c(0,200))+#,xlim=c(0,5))+
#                                             #scale_y_continuous(breaks=c(0,0.5,1.0))+
#                                             #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
#                                             theme_tufte()+
#                                             facet_grid(~protocol_fixed, labeller = label_parsed)+
#                                             my.theme#+
#                                             #theme(legend.title = element_text(colour="black", size=28, family="sans"),
#                                             	  #strip.text.y.right = element_text(angle = 0),
#                                             #	  axis.text.x=element_blank())
#                                             #legend.position = c(.95, .95),
#                                                   #legend.justification = c("right", "bottom"))#,
#                                                    # legend.margin = margin(6, 6, 6, 6))
                       
#                         ggsave(filename=paste0("tau_decay_boxplot"),plot=boxplot_tau_decay, device="jpeg",dpi=600, bg="white",units="in",width=plot_height,height=plot_width)

