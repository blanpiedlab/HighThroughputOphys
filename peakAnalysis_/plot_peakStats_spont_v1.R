### plot_peakStats_activity_vs_marker.R
## written by Samuel T Barlow 12.6.23

# Goal of this script is to accept a summary data.frame, peakStats, and plot data visualizations for Figure 2 of the optical physiology manuscript. 


subDir <- "spont_stats_v2"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/def_stim_vlines.R")
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_Figure2_v6.R")




plot_peakStats<- function(peak_stat_df = peak_stats, df = df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops){



  						#### establish variables that will be used for generating datasets





                        scalefactor = 0.75
                        cumdist_height = 16
                        cumdist_width = 3/5*cumdist_height
                        lineplot_width = 16
                        lineplot_height = 24
                        scatterplot_width = 16
                        scatterplot_height = 16
                        color_switch = !is.null(color_override)

                        interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)
                       



                        #### clean up dataframe by merging with peakStats

                        #peak_stat_df <- peak_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID) ) 
                        #ROIs_with_peaks <- unique(clean_df$ROI_key)


                        light_df <- df[,!(names(df) %in% remove_col)]
                        ROIs_to_remove <- light_df %>% group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%
                                                 dplyr::filter( is.na(dFF))
                        #print("These should be the borked ROIs we found (drift correction failures?)")
                        #print(unique(ROIs_to_remove$ROI_key))
                        #check<- ROIs_to_remove
                        ROIs_to_remove <- unique(ROIs_to_remove$ROI_key)                                                  


                        ROIs_with_stats <- unique(peak_stat_df$ROI_key)


                        clean_peak_stats<- peak_stat_df %>% dplyr::filter(!ROI_key %in% ROIs_to_remove) %>% mutate(peakgrp = paste0(ROI_key,"-",peakID))
                        clean_df_normTime <- suppressMessages(left_join(light_df, clean_peak_stats) %>% dplyr::filter(!ROI_key %in% ROIs_to_remove, ROI_key %in% ROIs_with_stats, peakID != "NotPeak") %>% 
                        					group_by_at(groupers) %>%
                        					mutate(peakgrp = paste0(ROI_key,"-",peakID),
                        						windowed_normTime = absoluteTime - absoluteTime[which.max(dFF)],
                        						max_windowTime = max(windowed_normTime),
                        						min_windowTime = min(windowed_normTime)) #%>%
                        					#dplyr::filter(max_windowTime < 1.0, min_windowTime > -1.0)
            							)
                        
                        print(unique(clean_df_normTime$peakgrp)[1:50])
                        getVariance <- suppressMessages(clean_df_normTime %>% group_by(Ca,peakgrp) %>%
                        					dplyr::filter(windowed_normTime > 0.04) %>%
                        					summarise(peak_variance = var(dFF,na.rm=TRUE))
                        					) 
                        print(getVariance)

                        peak_variance_scored <- getVariance %>% ungroup() #%>% dplyr::filter(peak_variance < quantile(peak_variance, 0.95, na.rm=TRUE))


                        clean_df_normTime <- clean_df_normTime %>% dplyr::filter(peakgrp %in% peak_variance_scored$peakgrp)
                        clean_peak_stats <- clean_peak_stats %>% dplyr::filter(peakgrp %in% peak_variance_scored$peakgrp)
     					
			new_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
			#segmentation_levels = c("marker","activity")

			clean_df_normTime$Ca<- as.factor(as.character( (clean_df_normTime$Ca) ) )
			clean_peak_stats$Ca<- as.factor(as.character( (clean_peak_stats$Ca) ) )
			

                        levels(clean_peak_stats$Ca) <- new_levels
                        levels(clean_df_normTime$Ca) <- new_levels

                    
                        


						#### generate ROI count from raw data

						getROIs<- clean_df_normTime %>% group_by(vid_key, segmentation,Ca, Ca_mM) %>% summarise(ROI_count = length(unique(ROI_key)))
						#getROIs$segmentation <- factor(getROIs$segmentation, levels = segmentation_levels)

						ROIs_with_peaks<- clean_df_normTime %>% group_by(segmentation,Ca, Ca_mM, ROI_key) %>% summarise(peaks_detected = length(unique(peakID))-1,
																		peak_positive = case_when(peaks_detected > 0 ~ 1,
																						peaks_detected == 0 ~ 0)
																		)


						area_imaging_region_um = 51.2 * 51.2
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





						scatter_ROIs<-ggplot(getROIs, aes(x=Ca_mM, y=ROI_count, group=segmentation, colour=segmentation))+#linetype=segmentation,
						                                            #geom_violin(size=1, alpha=0.2)+
						                                            stat_summary(geom="errorbar", fun.data="mean_se", size=1.5,alpha=1,width=(4/20) )+
						                                            stat_summary(geom="line", fun.y="mean", size=1.5,alpha=1)+
						                                            stat_summary(geom="point", fun.y="mean", size=5,alpha=1)+
						                                            labs( y="ROIs identified per video",
						                                                    x=expression(Ca^{"2+"}*" (mM)"),
						                                                    title="",#expression(Delta*"t (ms)"),	
						                                                    colour="Segmented by:",
					                                                    	tag = "F"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
						                                            #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
						                                            #coord_cartesian(ylim=c(0,200),xlim=c(0,6))+
						                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
						                                            coord_cartesian(xlim=c(0,4.1),ylim=c(0,30))+
						                                            scale_x_continuous(breaks = c(0.5,1,2,4))+
						                                            scale_y_continuous(breaks=c(0,10,20,30))+
						                            
						                                            theme_tufte()+
						                                            #facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(legend.title = element_text(colour="black", size=scalefactor*28, family="sans"),
						                                                  legend.position = c(0.5, 0.90),
					                                            	      strip.text.y.right = element_text(size=scalefactor*22, angle = 0),
					                                            	      axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans",angle=45, hjust=1))#,
						                                            	  #strip.text.y.right = element_text(angle = 0))
						                                            #legend.position = c(.95, .95),
						                                                  #legend.justification = c("right", "bottom"))#,
						                                                   # legend.margin = margin(6, 6, 6, 6))
						                       
							ggsave(filename=paste0("scatter_ROIs.jpeg"),plot=scatter_ROIs, device="jpeg",dpi=600, bg="white",units="in",width=10,height=10)

							ROI_per_vid  <<- scatter_ROIs

							scatter_ROIs_with_peaks<-ggplot(ROIs_with_peaks, aes(x=Ca_mM, y=peak_positive, group=segmentation, colour=segmentation))+#linetype=segmentation,
						                                            #geom_violin(size=1, alpha=0.2)+
						                                            stat_summary(geom="errorbar", fun.data="mean_se", size=1.5,alpha=1,width=(4/20) )+
						                                            stat_summary(geom="line", fun.y="mean", size=1.5,alpha=1)+
						                                            stat_summary(geom="point", fun.y="mean", size=5,alpha=1)+
						                                            labs( y="Probability of iGluSnFR3 activity",
						                                                    x=expression(Ca^{"2+"}*" (mM)"),
						                                                    title="",#expression(Delta*"t (ms)"),	
						                                                    #colour="Segmented by:",
					                                                    	tag = "G"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
						                                            coord_cartesian(xlim=c(0,4.1),ylim=c(0,1))+
						                                            scale_x_continuous(breaks = c(0.5,1,2,4))+
						                                            scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
						                            
						                                            theme_tufte()+
						                                            #facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(#legend.title = element_text(colour="black", size=28, family="sans"),
						                                                  legend.position = "none",
					                                            	      #strip.text.y.right = element_text(size=22, angle = 0),
					                                            	      axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans",angle=45, hjust=1))#,
						                                            	  #strip.text.y.right = element_text(angle = 0))
						                                            #legend.position = c(.95, .95),
						                                                  #legend.justification = c("right", "bottom"))#,
						                                                   # legend.margin = margin(6, 6, 6, 6))
						                       
							ggsave(filename=paste0("scatter_ROIs_with_peaks.jpeg"),plot=scatter_ROIs_with_peaks, device="jpeg",dpi=600, bg="white",units="in",width=10,height=10)

							peak_probability  <<- scatter_ROIs_with_peaks







 ####### plot cumdists
                         xmax = max(clean_peak_stats$amplitude)
                        ecdf_amplitude<-ggplot(clean_peak_stats, aes(x=amplitude, group=Ca, linetype=segmentation,colour=Ca))+
                                             stat_ecdf(size=2,alpha=1, lineend="round")+
                                             
                                            labs( x=expression("Peak "*Delta*"F/F"),
                                                    y="Cumulative Density",
                                                    title="Peak amplitude"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override, guide= "none")}+
                                            scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,1.1),xlim=c(0,2.5))+
                                            scale_y_continuous(breaks=c(0,0.5,1.0))+
                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            #facet_grid(Ca~., labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.position = "none",
                                            	  strip.text.y.right = element_text(angle = 0))
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("amplitude_cumdist.jpeg"),plot=ecdf_amplitude, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)



                        xmax = max(clean_peak_stats$tau_decay_ms)
                         
                         ecdf_tau<-ggplot(clean_peak_stats, aes(x=tau_decay_ms, group=Ca, linetype=segmentation,colour=Ca))+
                                             stat_ecdf(size=2,alpha=1, lineend="round")+
                                             
                                            labs( x=expression(tau[decay]*" (ms)"),
                                                    y="Cumulative Density",
                                                    title=expression(tau[decay])
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override, guide= "none")}+
                                            scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,1.1),xlim=c(0,200))+
                                            scale_y_continuous(breaks=c(0,0.5,1.0))+
                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            #facet_grid(Ca~., labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.position = "none",
                                            	  strip.text.y.right = element_text(angle = 0))
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("tau_decay_ms_cumdist.jpeg"),plot=ecdf_tau, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)



                        #xmax = max(clean_peak_stats$interSpike_ms)
                         

                        #  ks_test<- wilcox.test(clean_peak_stats$amplitude[which(clean_peak_stats$segmentation == "marker",clean_peak_stats$Ca == "1Ca")],clean_peak_stats$amplitude[which(clean_peak_stats$segmentation == "activity",clean_peak_stats$Ca == "1Ca")])
                        #  print("testing amplitude with wilcox.test")
                        #  print(ks_test) 
                         
                        #  ks_test<- wilcox.test(clean_peak_stats$tau_decay_ms[which(clean_peak_stats$segmentation == "marker",clean_peak_stats$Ca == "1Ca")],clean_peak_stats$tau_decay_ms[which(clean_peak_stats$segmentation == "activity",clean_peak_stats$Ca == "1Ca")])
                        #  print("testing tau_decay_ms with wilcox.test")
                        #  print(ks_test) 
                         
                        


                        # ks_test<- wilcox.test(clean_peak_stats$t_rise[which(clean_peak_stats$segmentation == "marker",clean_peak_stats$Ca == "1Ca")],clean_peak_stats$t_rise[which(clean_peak_stats$segmentation == "activity",clean_peak_stats$Ca == "1Ca")])
                        #  print("testing t_rise with wilcox.test")
                        #  print(ks_test) 
                         
                        #  ks_test<- wilcox.test(clean_peak_stats$t_decay[which(clean_peak_stats$segmentation == "marker",clean_peak_stats$Ca == "1Ca")],clean_peak_stats$t_decay[which(clean_peak_stats$segmentation == "activity",clean_peak_stats$Ca == "1Ca")])
                        #  print("testing t_decay with wilcox.test")
                        #  print(ks_test) 
                         
                        #  ks_test<- wilcox.test(clean_peak_stats$t_half[which(clean_peak_stats$segmentation == "marker",clean_peak_stats$Ca == "1Ca")],clean_peak_stats$t_half[which(clean_peak_stats$segmentation == "activity",clean_peak_stats$Ca == "1Ca")])
                        #  print("testing t_half with wilcox.test")
                        #  print(ks_test) 
                         



### as boxplots
    				violin_amplitude<-ggplot(clean_peak_stats, aes(x=segmentation, y=amplitude, group=segmentation, colour=Ca))+#linetype=segmentation,
                                             #stat_ecdf(size=2,alpha=1, lineend="round")+
                                             geom_sina(alpha=0.3)+
                                             geom_violin(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75))+
                                            labs( y=expression("Peak "*Delta*"F/F"),
                                                    x="",#"Cumulative Density",
                                                    title="Peak amplitude",
                                                    colour="Segmented by:"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,2.5))+#,xlim=c(0,5))+
                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            facet_grid(~Ca, labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.title = element_text(colour="black", size=28, family="sans"),
                                            	  #strip.text.y.right = element_text(angle = 0),
                                            	  axis.text.x=element_blank())
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("amplitude_violin.jpeg"),plot=violin_amplitude, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)

                        violin_tau_decay_ms<-ggplot(clean_peak_stats, aes(x=segmentation, y=tau_decay_ms, group=segmentation, colour=Ca))+#linetype=segmentation,
                                             #stat_ecdf(size=2,alpha=1, lineend="round")+
                                             geom_sina(alpha=0.3)+
                                             geom_violin(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75))+
                                            labs( y=expression(tau[decay]*" (ms)"),
                                                    x="",#"Cumulative Density",
                                                    title=expression(tau[decay]*" (ms)"),
                                                    colour="Segmented by:"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,200))+#,xlim=c(0,5))+
                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            facet_grid(~Ca, labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.title = element_text(colour="black", size=28, family="sans"),
                                            	  #strip.text.y.right = element_text(angle = 0),
                                            	  axis.text.x=element_blank())
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("tau_decay_ms_violin.jpeg"),plot=violin_tau_decay_ms, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)



                  

##### HISTOGRAMS
						  lower_bound_x = 0
	                      upper_bound_x = 10
      
                         bin = (upper_bound_x - lower_bound_x) / 50  


                        cleanup_time_data <- clean_peak_stats %>% dplyr::filter(t_rise > 1 & !is.na(t_rise)) %>% 
                           										dplyr::filter(t_half > 1 & !is.na(t_half)) %>%
                           										dplyr::filter(t_decay > 1 & !is.na(t_decay))
                        # time_data_medians<- cleanup_time_data %>%
                        # 					group_by(Ca, segmentation) %>%
                        # 					summarise(median_t_rise = median(t_rise,na.rm=TRUE),
			# 								  			median_t_half = median(t_half,na.rm=TRUE),
			# 								  			median_t_decay = median(t_decay,na.rm=TRUE) )
		                           



                        my_medians <- clean_peak_stats %>%
									  	group_by(Ca, segmentation) %>%
									  summarize(median_amplitude = median(amplitude),
									  			median_tau = median(tau_decay_ms),
									  			#median_interSpike =median(interSpike_ms),
									  			median_t_rise = median(t_rise,na.rm=TRUE),
											  			median_t_half = median(t_half,na.rm=TRUE),
											  			median_t_decay = median(t_decay,na.rm=TRUE) ) #,

					   histo_amplitude<-ggplot(clean_peak_stats, aes(x=amplitude, group=segmentation, fill=Ca,colour=Ca))+
					                                            geom_vline(data=my_medians, aes(xintercept = median_amplitude, colour=Ca), size=scalefactor*2, lty="solid" )+
                                        						
					                                            geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression("Peak "*Delta*"F/F"),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "I"
					                                                    )+
					                                            
					                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                		    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(labels=new_levels, values=color_override)}+
                                                                		    #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
					                                            coord_cartesian(ylim=c(0,1.3),xlim=c(0,upper_bound_x))+
					                                            scale_y_continuous(breaks=c(0,0.5,1.0), expand=c(0,0))+
					                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
					                            
					                                            theme_tufte()+
					                                            facet_grid(Ca~., labeller = label_parsed)+
					                                            my.theme+
					                                            theme(legend.position = "none",
					                                            	  legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
					                                            	  strip.text.y.right = element_blank())
					                                            #legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_amplitude.jpeg"),plot=histo_amplitude, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)

					     #amplitude_histo <<- histo_amplitude




					       lower_bound_x = 0
	                      upper_bound_x = 200
      
                         bin = (upper_bound_x - lower_bound_x) / 50  

					     histo_tau_decay_ms<-ggplot(clean_peak_stats, aes(x=tau_decay_ms, group=segmentation, fill=Ca,colour=Ca))+
					                                            geom_vline(data=my_medians, aes(xintercept = median_tau, colour=Ca), size=scalefactor*2, lty="solid" )+
                                        						
					                                            geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					 
					                                            labs( x=expression(tau[decay]*" (ms)"),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#""expression(tau[decay]*" (ms)"),
					                                                    fill="Segmented by:",
					                                                    tag = "J"
					                                                    )+
					                                            
					                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                		    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(labels=new_levels, values=color_override)}+
                                                                		    #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
					                                            coord_cartesian(ylim=c(0,1.3),xlim=c(0,upper_bound_x))+
					                                            scale_y_continuous(breaks=c(0,0.5,1.0), expand=c(0,0))+
					                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
					                            
					                                            theme_tufte()+
					                                            facet_grid(Ca~., labeller = label_parsed)+
					                                            my.theme+
					                                            
					                                            theme(legend.position = "none",
					                                            	  legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
					                                            	  strip.text.y.right = element_blank())#legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_tau_decay_ms.jpeg"),plot=histo_tau_decay_ms, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)

					     #tau_histo <<- histo_tau_decay_ms

					       lower_bound_x = 0
	                      upper_bound_x = 150
      
                         bin = (upper_bound_x - lower_bound_x) / 50  

					      



					       lower_bound_x = 0
	                      upper_bound_x = 150
      
                         bin = (upper_bound_x - lower_bound_x) / 50  



					     histo_t_rise<-ggplot(cleanup_time_data, aes(x=t_rise, group=segmentation, fill=Ca,colour=Ca))+
					                                            geom_vline(data=my_medians, aes(xintercept = median_t_rise, colour=Ca), size=scalefactor*2, lty="solid" )+
                                        						
					                                            geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression(t[rise]),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "I"
					                                                    )+
					                                            
					                                           {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                		    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(labels=new_levels, values=color_override)}+
                                                                		    #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
					                                            coord_cartesian(ylim=c(0,1.3),xlim=c(0,upper_bound_x))+
					                                            scale_y_continuous(breaks=c(0,0.5,1.0), expand=c(0,0))+
					                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
					                            
					                                            theme_tufte()+
					                                            facet_grid(Ca~., labeller = label_parsed)+
					                                            my.theme+
					                                            theme(legend.position = "none",
					                                            	  legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
					                                            	  strip.text.y.right = element_blank())
					                                            #legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_t_rise.jpeg"),plot=histo_t_rise, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)


					        lower_bound_x = 0
	                      upper_bound_x = 400
      
                         bin = (upper_bound_x - lower_bound_x) / 50 

					     histo_t_decay<-ggplot(cleanup_time_data, aes(x=t_decay, group=segmentation, fill=Ca,colour=Ca))+
					                                            geom_vline(data=my_medians, aes(xintercept = median_t_decay, colour=Ca), size=scalefactor*2, lty="solid" )+
                                        						
					                                            geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression(t[decay]),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "I"
					                                                    )+
					                                            
					                                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                		    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(labels=new_levels, values=color_override)}+
                                                                		    #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
					                                            coord_cartesian(ylim=c(0,1.3),xlim=c(0,upper_bound_x))+
					                                            scale_y_continuous(breaks=c(0,0.5,1.0), expand=c(0,0))+
					                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
					                            
					                                            theme_tufte()+
					                                            facet_grid(Ca~., labeller = label_parsed)+
					                                            my.theme+
					                                            theme(legend.position = "none",
					                                            	  legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
					                                            	  strip.text.y.right = element_blank())
					                                            #legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_t_decay.jpeg"),plot=histo_t_decay, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)


					        lower_bound_x = 0
	                      upper_bound_x = 200
      
                         bin = (upper_bound_x - lower_bound_x) / 50 

					     histo_t_half<-ggplot(cleanup_time_data, aes(x=t_half, group=segmentation, fill=Ca,colour=Ca))+
					                                            geom_vline(data=my_medians, aes(xintercept = median_t_half, colour=Ca), size=scalefactor*2, lty="solid" )+
                                        						
					                                            geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression(t[1/2]),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "I"
					                                                    )+
					                                            
					                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                		    {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(labels=new_levels, values=color_override)}+
                                                                		    #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            coord_cartesian(ylim=c(0,1.3),xlim=c(0,upper_bound_x))+
					                                            scale_y_continuous(breaks=c(0,0.5,1.0), expand=c(0,0))+
					                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
					                            
					                                            theme_tufte()+
					                                            facet_grid(Ca~., labeller = label_parsed)+
					                                            my.theme+
					                                            theme(legend.position = "none",
					                                            	  legend.title = element_blank(),#element_text(colour="black", size=28, family="sans"),
					                                            	  strip.text.y.right = element_blank())
					                                            #legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_t_half.jpeg"),plot=histo_t_half, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)


####### SCATTERPLOTS ####### 




                        plot_avgs <- clean_peak_stats %>%
									  	group_by(Ca, segmentation) %>%
									  summarize(n = n(),
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

									        # mean_dt = mean(interSpike_ms),
									         #sd_dt = sd(interSpike_ms),
									         #se_dt = sd_dt / sqrt(n),
									         )#,
									         
									         #ci = qt(0.975, df = n - 1) * sd / sqrt(n))

					   print(plot_avgs)

					   	abridged_peak_stats <- clean_peak_stats %>% dplyr::filter(tau_decay_ms <=75) %>% group_by(Ca) %>% 
					   												mutate(ampGrp_vesicles = case_when(amplitude < 0.7 ~ 1,
					   																			amplitude >= 0.7 & amplitude < 1.5 ~ 2,
					   																			amplitude >= 1.5 & amplitude < 2.4 ~ 3,
					   																			amplitude >= 2.4 & amplitude < 3.8 ~ 4,
					   																			amplitude >= 3.8 & amplitude < 9 ~ 5)
					   												)

					   	abridged_bin_one<- abridged_peak_stats %>% dplyr::filter(ampGrp_vesicles == 1)
						abridged_bin_two<- abridged_peak_stats %>% dplyr::filter(ampGrp_vesicles == 2)

					   	total_rows_bin1<- nrow(abridged_bin_one)
					   	total_rows_bin2<- nrow(abridged_bin_two)
					   	

					   	rand_in <- sample(total_rows_bin1,total_rows_bin1/10,replace=FALSE)
					   	resampled_bin_one<- abridged_bin_one[rand_in,]
					   	# rand_in <- sample(total_rows_bin2,total_rows_bin2/10,replace=FALSE)
					   	
					   	# resampled_bin_two<-abridged_bin_two[rand_in,]
					   	abridged_other_bins<- abridged_peak_stats %>% dplyr::filter(ampGrp_vesicles != 1)#, ampGrp_vesicles != 2)

					   	resampled_peak_stats<- bind_rows(resampled_bin_one,abridged_other_bins)
					   	check_peak_stats<- abridged_peak_stats %>% group_by(Ca, ampGrp_vesicles) %>% tally() 
					   	as_percentage <- check_peak_stats %>% group_by(Ca) %>% summarise(ampGrp = ampGrp_vesicles,
					   																		n = n,
					   																		n_total= sum(n),
					   																	percentage = n/n_total*100)
					   	print(check_peak_stats)

						scatter_amplitude_tau<-ggplot(resampled_peak_stats, aes(x=amplitude/0.4, y=tau_decay_ms, group=segmentation, colour=Ca))+#linetype=segmentation,
						                                            geom_vline(xintercept=c(1,2,5,10,20),lty="dashed",colour="black")+
						                                            
						                                            geom_point(size=2, alpha=0.5)+
						                                            #geom_errorbar(data=plot_avgs, aes(x=mean_amplitude,y=mean_tau, ymin = mean_tau - se_tau, ymax = mean_tau + se_tau, group=segmentation, colour=segmentation, width=0.15), alpha=1)+
						                                            #geom_errorbarh(data=plot_avgs, aes(x=mean_amplitude, y=mean_tau, xmin = mean_amplitude - se_tau, xmax = mean_amplitude + se_tau, colour=segmentation,height=5), alpha=1)+
						                                            
						                                            #geom_point(data=plot_avgs, aes(x=mean_amplitude, y=mean_tau, group=segmentation, colour=segmentation),size=2, alpha=1)+
						                                            labs( y=expression(tau[decay]*" (ms)"),
						                                                    x=expression(Evoked[Delta*"F/F"] / Spont[Delta*"F/F"]),#expression("Peak "*Delta*"F/F"),
						                                                    title="",#expression(Delta*"t (ms)"),
						                                                    colour="Segmented by:"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
                                                                #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
						                                            coord_cartesian(ylim=c(0,100),xlim=c(0,(8/0.4)))+
						                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
						                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
						                            
						                                            theme_tufte()+
						                                            facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(legend.title = element_text(colour="black", size=28, family="sans"),
						                                                  legend.position = "none",
					                                            	      strip.text.y.right = element_text(size=22, angle = 0))#,
						                                            	  #strip.text.y.right = element_text(angle = 0))
						                                            #legend.position = c(.95, .95),
						                                                  #legend.justification = c("right", "bottom"))#,
						                                                   # legend.margin = margin(6, 6, 6, 6))
						                       
							ggsave(filename=paste0("scatter_amplitude_tau_as_ratio.jpeg"),plot=scatter_amplitude_tau, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width*1.5,height=cumdist_height)






										   #summarise(overall_peak_freq = mean(peak_freq,na.rm=TRUE))




##### AVG PPR WAVEFORMS ######


						avgtrace_Plot<-ggplot(clean_df_normTime, aes(x=windowed_normTime, y = dFF,  group=peakgrp, colour=Ca, linetype = segmentation))+
						                            geom_path(colour="grey61",size=0.3,alpha=0.3)+
						                            stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.8)+
						                            
						                           	labs( x="Normalized time (s)",
						                                   y=expression(Delta*"F/F") )+#,
						                                  #colour="Segmented by:")+
						                                  # #colour="Segmented by:",
					                                          #           tag = "H")+
						                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
									    scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						                            scale_x_continuous(breaks = c(0,0.1,0.2))+#,limits=c(-0.05,0.2))+
						                            scale_y_continuous(breaks=c(0,1.0,2,3,4,5,6,7))+
						                            theme_tufte()+
						                            #guides(colour="none",fill="none")+
						                            #coord_cartesian(ylim=c(-0.1,1.0),xlim=c(-0.05,0.2))+
									    my.theme+
						                            facet_grid(~Ca, labeller = label_parsed)#+
						                            #theme(legend.title = element_text(colour="black", size=28, family="sans"),
						                           # 	legend.position="none")

						                            #avg_peaks  <<- avg_PP_tracePlot
						                             #       rm(avg_PP_tracePlot
						                                

						                             ggsave(filename=paste0("averagePeaks_.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)
                         

						#avg_traces <<- avgtrace_Plot
					     
						#long_peaks<- clean_peak_stats %>% dplyr::filter(tau_decay_ms > 200)
						clean_peaks_long<- clean_df_normTime %>%  dplyr::filter(amplitude >= 1.2, Ca == "\"4 mM \" * Ca^\"2+\"")#t_decay > 50, t_decay < 150)dplyr::filter(tau_decay_ms > 50)




						 avgtrace_Plot<-ggplot(clean_peaks_long, aes(x=windowed_normTime, y = dFF,  group=peakgrp, colour=Ca, linetype = segmentation))+
						                            geom_path(colour="grey61",size=0.5,alpha=0.3)+
						                            stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.8)+
						                            geom_hline(yintercept = 0.4, lty="dashed",colour="black")+
						                           	labs( x="Normalized time (s)",
						                                  y=expression(Delta*"F/F") )+#,
						                                  #colour="Segmented by:")+
						                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
									    scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						                            scale_x_continuous(breaks = c(-0.25,0,0.25))+#,limits=c(-0.05,0.2))+
						                            scale_y_continuous(breaks=c(0,0.5,1.0,2,3,4,5,6,7))+
						                            theme_tufte()+
						                            #guides(colour="none",fill="none")+
						                           #coord_cartesian(ylim=c(-0.1,1.0),xlim=c(-0.05,0.2))+
									    my.theme+
						                            facet_grid(~Ca, labeller = label_parsed)#+
						                            #theme(legend.title = element_text(colour="black", size=28, family="sans"))

						                            #avg_peaks  <<- avg_PP_tracePlot
						                             #       rm(avg_PP_tracePlot)
						                                

						                             ggsave(filename=paste0("averagePeaks_long_cutoff_above_4dFF.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)

						   



						    clean_peaks_short<- clean_df_normTime %>% dplyr::filter(amplitude < 2.4 & amplitude >=1.2)#dplyr::filter(t_decay < 50)  dplyr::filter(tau_decay_ms <= 50 )#





						   avgtrace_Plot<-ggplot(clean_peaks_short, aes(x=windowed_normTime, y = dFF,  group=peakgrp, colour=Ca, linetype = segmentation))+
						                            geom_path(colour="grey61",size=0.5,alpha=0.3)+
						                            stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.8)+
						                            
						                           	labs( x="Normalized time (s)",
						                                   y=expression(Delta*"F/F") )+#,
						                                  #colour="Segmented by:")+
						                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override)}+
									    scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						                            scale_x_continuous(breaks = c(-0.25,0,0.25))+#,limits=c(-0.05,0.2))+
						                            scale_y_continuous(breaks=c(0,0.5,1.0,2,3,4,5,6,7))+
						                            theme_tufte()+
						                            #guides(colour="none",fill="none")+
						                            #coord_cartesian(ylim=c(-0.1,1.0),xlim=c(-0.05,0.2))+
									    my.theme+
						                            facet_grid(~Ca, labeller = label_parsed)#+
						                            #theme(legend.title = element_text(colour="black", size=28, family="sans"))

						                            #avg_peaks  <<- avg_PP_tracePlot
						                             #       rm(avg_PP_tracePlot)
						                                

						                             ggsave(filename=paste0("averagePeaks_short_cutoff_bin3.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)

list.output = list(ROI_totals = getROIs, ROI_probability = ROIs_with_peaks, ROIs_per_vid_summary = check_ROIs_per_vid, Medians = my_medians, Averages = plot_avgs, clean_stats =clean_peak_stats,clean_df = clean_df_normTime,vesicle_estimate = check_peak_stats,vesicle_estimate_pct = as_percentage)
list.output
}





