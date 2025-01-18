### plot_peakStats_activity_vs_marker.R
## written by Samuel T Barlow 12.6.23

# Goal of this script is to accept a summary data.frame, peakStats, and plot data visualizations for Figure 2 of the optical physiology manuscript. 


subDir <- "peakStats_activity_vs_marker_v3"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/def_stim_vlines.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_Figure2_v6.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_v4.R")




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


                        

                        peak_stat_df <- peak_stat_df %>% mutate(windowedPeakID = paste0("windowed",peakID),
                        				 imaging_region = gsub("-repl\\d\\d", "", vid_key)) 

                        light_df <- df[,!(names(df) %in% remove_col)] %>% 
                        				 mutate(imaging_region = gsub("-repl\\d\\d", "", vid_key)) 


                        
                         ROIs_to_remove_timeskip<- light_df %>% #only retain ROIs that are do not have a frame skip in time of interest
                                                 group_by(vid_key, Ca, Ca_mM,trackROI_key,ROI_key, ROINumber, exposeNum) %>%            
                                                 summarise(maxInterframe = max(interFrame, na.rm=TRUE)) %>%
                                                 dplyr::filter(maxInterframe > 0.18) %>%
                                                 select(ROI_key,maxInterframe) 
                        ROIs_to_remove_timeskip<-unique(ROIs_to_remove_timeskip$ROI_key)                                             



                        clean_df <- light_df %>%  dplyr::filter( !ROI_key %in% ROIs_to_remove_timeskip) 

                        clean_peaks<- suppressMessages(left_join(clean_df,peak_stat_df) %>% dplyr::filter(windowedPeakID != "NotPeak") ) 
 
                        clean_peaks_normTime<- suppressMessages(clean_peaks %>% group_by_at(groupers) %>% mutate(time_max = absoluteTime[which.max(dFF)],
                        										      new_normTime = absoluteTime - time_max) ) 

     					
						new_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
						segmentation_levels = c("activity", "marker")
                        
     					clean_peaks_normTime$Ca<- as.factor(as.character( (clean_peaks_normTime$Ca) ) )
     					peak_stat_df$Ca<- as.factor(as.character( (peak_stat_df$Ca) ) )
     					

                        levels(peak_stat_df$Ca) <- new_levels
                        levels(clean_peaks_normTime$Ca) <- new_levels
                        clean_peaks_normTime$segmentation <- factor(clean_peaks_normTime$segmentation, levels = segmentation_levels)


                    
                        


						#### generate ROI count from raw data

						getROIs<- clean_df %>% #dplyr::filter(!ROI_key %in% ROIs_to_remove) %>% 
									group_by(vid_key, segmentation,Ca, Ca_mM) %>% summarise(ROI_count = length(unique(ROI_key)))
						getROIs$segmentation <- factor(getROIs$segmentation, levels = rev(segmentation_levels))

						ROIs_with_peaks<- clean_df %>% 
										#dplyr::filter(!ROI_key %in% ROIs_to_remove,stimEpoch == 1) %>% 
										group_by(segmentation,Ca, Ca_mM, ROI_key) %>% summarise(peaks_detected = length(unique(peakID))-1,
																									peak_positive = case_when(peaks_detected > 0 ~ 1,
																									peaks_detected == 0 ~ 0)
																									)
						ROIs_with_peaks$segmentation <- factor(ROIs_with_peaks$segmentation, levels = segmentation_levels)
						


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




#### Generate 4 raster plots of different ROIs. 

###### for loop for plotting spatial data
        png_dir = "Y:\\Sam/paper1_datasets/spont_v1/imageFiles/zstack/zstackOutputs_SQ"#/GluSnFR3_SyPhy_dish06_plate06_region04_2Ca_zstack_lowpwr__Results"
        prefix_str = "^MAX_GluSnFR3_SyPhy_dish06_plate06_region04_2Ca"
        suffix_str = "_zstack_lowpwr__FLATTENED.png"
        all_lut_str = "_zstack_lowpwr__FLATTENED_COLOR.png"
        dir_regex = "GluSnFR3_SyPhy_dish06_plate06_region04"

  		print("Taking a crack at generating png plots")

        subset_df<- clean_df %>% mutate(formatted_vid_key = str_extract(vid_key, "dish\\d\\d-plate\\d\\d-region\\d\\d")) %>% dplyr::filter(formatted_vid_key %in% c("dish05-plate01-region02"))

        
        vid_key_vec <- unique(subset_df$formatted_vid_key)
        print(vid_key_vec)
        save_bool=TRUE
        #guide_bool = TRUE
        guide_title = expression("Segmented by:") # expression(CV[tau*"decay"])#
        
        tmp_tag_var_override = c("B","B","B","B")

        marker_map = "marker_ROImap.zip"

        ###scale bar params! set the left-corner 
        user_x = 20
        user_y = 24.5




        for (k in 1:length(vid_key_vec)){
            
                vid_df<- subset_df %>% mutate(fix_Ca = paste0(Ca_mM, "Ca") )
                Ca_vec <- unique(vid_df$fix_Ca)
                print(Ca_vec)
            	special_tag = "A"    
                tmp_plot = png_plotter_v4(df = vid_df, png_dir = png_dir,prefix_str=prefix_str, suffix_str = suffix_str,
                                        all_lut_str = all_lut_str, dir_regex = dir_regex, vars_to_plot=vars_to_plot,save_bool=save_bool,
                                        tag_var = special_tag,user_x = user_x, user_y = user_y)
                   nam <- paste("map_color_", special_tag, sep = "")
		   assign(nam, tmp_plot,envir = .GlobalEnv)

                count=0
                
                activity_map = "activity_ROImap.zip"
                save_bool = TRUE
                guide_bool = TRUE
                tmp_tag_var = "B"
                i=2
                #                     
                # for(i in 1:length(Ca_vec)) {
                # 					current_Ca = Ca_vec[i]
                # 					print(paste0("Current_Ca is == ",current_Ca))
                # 					#if(current_Ca == "0.5Ca") {activity_map = "activity_ROImap.zip"} 
                # 					#if(current_Ca == "1Ca") { activity_map = "activity_ROImap.zip"}
                # 					if(current_Ca == "2Ca") { activity_map = "activity_ROImap.zip"}
                # 					#if(current_Ca == "4Ca") { activity_map = "activity_ROImap.zip"}
                #                     tmp_tag_var = tmp_tag_var_override[i]
                #                     tmp_plot = NULL
                #                     Ca_df<- vid_df %>% dplyr::filter(fix_Ca == current_Ca)
                #                     if(save_bool==TRUE) { 
                                        

                #                        if(i==2){guide_bool = TRUE} else {guide_bool = FALSE}

                                       tmp_plot = png_plotter_ROIoverlay(df = vid_df, 
                                       									png_dir = png_dir,
                                       									prefix_str=prefix_str, 
                                       									suffix_str = suffix_str, 
                                       									dir_regex = dir_regex, 
                                       									marker_map = marker_map, 
                                       									activity_map=activity_map, 
                                       									save_bool=save_bool,
                                       									tag_var=tmp_tag_var, 
                                       									user_x = user_x, 
                                       									user_y = user_y,
                                       									guide_title=guide_title,
                                       									guide_bool=guide_bool,
                                       									scalefactor=scalefactor)

                                       	nam <- paste("map_", tmp_tag_var_override[i], sep = "")
                                        assign(nam, tmp_plot,envir = .GlobalEnv)
#                                    }                                

                #}

        }
                                    























						scatter_ROIs<-ggplot(getROIs, aes(x=Ca_mM, y=ROI_count, group=segmentation, colour=segmentation))+#linetype=segmentation,
						                                            #geom_violin(size=1, alpha=0.2)+
						                                            stat_summary(geom="errorbar", fun.data="mean_se", size=1,alpha=1,width=(4/20) )+
						                                            stat_summary(geom="line", fun.y="mean", size=1,alpha=1)+
						                                            stat_summary(geom="point", fun.y="mean", size=3,alpha=1)+
						                                            labs( y="ROIs per video",
						                                                    x=expression(Ca^{"2+"}*" (mM)"),
						                                                    title="",#expression(Delta*"t (ms)"),	
						                                                    #colour="Segmented by:",
					                                                    	tag = "C"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
						                                            #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
						                                            #coord_cartesian(ylim=c(0,200),xlim=c(0,6))+
						                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
						                                            coord_cartesian(xlim=c(0,4.1),ylim=c(0,30))+
						                                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
						                                            scale_y_continuous(breaks=c(0,10,20,30))+
						                            
						                                            theme_tufte()+
						                                            #facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(legend.title = element_blank(),#element_text(colour="black", size=scalefactor*34, family="sans"),
						                                                  legend.position = "top")#,
					                                            	      #strip.text.y.right = element_text(size=scalefactor*22, angle = 0),
					                                            	      #axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans",angle=45, hjust=1))#,
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
					                                                    	tag = "E"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
						                                            coord_cartesian(xlim=c(0,4.1),ylim=c(0,1))+
						                                            scale_x_continuous(breaks=c(0.5,1,2,4),labels=c("0.5","1","2","4"))+
                            								    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
						                            
						                                            theme_tufte()+
						                                            #facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(#legend.title = element_text(colour="black", size=28, family="sans"),
						                                                  legend.position = "none")#,
					                                            	      #strip.text.y.right = element_text(size=22, angle = 0),
					                                            	      #axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans",angle=45, hjust=1))#,
						                                            	  #strip.text.y.right = element_text(angle = 0))
						                                            #legend.position = c(.95, .95),
						                                                  #legend.justification = c("right", "bottom"))#,
						                                                   # legend.margin = margin(6, 6, 6, 6))
						                       
							ggsave(filename=paste0("scatter_ROIs_with_peaks.jpeg"),plot=scatter_ROIs_with_peaks, device="jpeg",dpi=600, bg="white",units="in",width=10,height=10)

							peak_probability  <<- scatter_ROIs_with_peaks







 ####### plot cumdists
                         xmax = max(peak_stat_df$amplitude)
                        ecdf_amplitude<-ggplot(peak_stat_df, aes(x=amplitude, group=segmentation, linetype=segmentation,colour=Ca))+
                                             stat_ecdf(size=2,alpha=1, lineend="round")+
                                             
                                            labs( x=expression("Peak "*Delta*"F/F"),
                                                    y="Cumulative Density",
                                                    title="Peak amplitude"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=new_levels, values=color_override, guide= "none")}+
                                            scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,1.1),xlim=c(0,6))+
                                            scale_y_continuous(breaks=c(0,0.5,1.0))+
                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            facet_grid(Ca~., labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.position = "none",
                                            	  strip.text.y.right = element_text(angle = 0))
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("amplitude_cumdist.jpeg"),plot=ecdf_amplitude, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)



                        xmax = max(peak_stat_df$tau_decay_ms)
                         
                         ecdf_tau<-ggplot(peak_stat_df, aes(x=tau_decay_ms, group=segmentation, linetype=segmentation,colour=Ca))+
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
                                            facet_grid(Ca~., labeller = label_parsed)+
                                            my.theme+
                                            theme(legend.position = "none",
                                            	  strip.text.y.right = element_text(angle = 0))
                                            #legend.position = c(.95, .95),
                                                  #legend.justification = c("right", "bottom"))#,
                                                   # legend.margin = margin(6, 6, 6, 6))
                       
                        ggsave(filename=paste0("tau_decay_ms_cumdist.jpeg"),plot=ecdf_tau, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)




                        #   ks_test<- wilcox.test(peak_stat_df$amplitude[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$amplitude[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        #   print("testing amplitude with wilcox.test")
                        #   print(ks_test) 
                         
                        #   ks_test<- wilcox.test(peak_stat_df$tau_decay_ms[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$tau_decay_ms[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        #   print("testing tau_decay_ms with wilcox.test")
                        #   print(ks_test) 
                         
                        # #  ks_test<- wilcox.test(peak_stat_df$interSpike_ms[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$interSpike_ms[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        # #  print("testing interSpike_ms with wilcox.test")
                        # #  print(ks_test) 
                        


                        #  ks_test<- wilcox.test(peak_stat_df$t_rise[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$t_rise[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        #   print("testing t_rise with wilcox.test")
                        #   print(ks_test) 
                         
                        #   ks_test<- wilcox.test(peak_stat_df$t_decay[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$t_decay[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        #   print("testing t_decay with wilcox.test")
                        #   print(ks_test) 
                         
                        #   ks_test<- wilcox.test(peak_stat_df$t_half[which(peak_stat_df$segmentation == "marker",peak_stat_df$Ca == "1Ca")],peak_stat_df$t_half[which(peak_stat_df$segmentation == "activity",peak_stat_df$Ca == "1Ca")])
                        #   print("testing t_half with wilcox.test")
                        #   print(ks_test) 
                         



### as boxplots
    				violin_amplitude<-ggplot(peak_stat_df, aes(x=segmentation, y=amplitude, group=segmentation, colour=segmentation))+#linetype=segmentation,
                                             #stat_ecdf(size=2,alpha=1, lineend="round")+
                                             geom_sina(alpha=0.3)+
                                             geom_violin(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75))+
                                            labs( y=expression("Peak "*Delta*"F/F"),
                                                    x="",#"Cumulative Density",
                                                    title="Peak amplitude",
                                                    colour="Segmented by:"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
                                            #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
                                            coord_cartesian(ylim=c(0,6))+#,xlim=c(0,5))+
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

                        violin_tau_decay_ms<-ggplot(peak_stat_df, aes(x=segmentation, y=tau_decay_ms, group=segmentation, colour=segmentation))+#linetype=segmentation,
                                             #stat_ecdf(size=2,alpha=1, lineend="round")+
                                             geom_sina(alpha=0.3)+
                                             geom_violin(fill=NA,draw_quantiles = c(0.25, 0.5, 0.75))+
                                            labs( y=expression(tau[decay]*" (ms)"),
                                                    x="",#"Cumulative Density",
                                                    title=expression(tau[decay]*" (ms)"),
                                                    colour="Segmented by:"
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
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
	                      upper_bound_x = 2
      
                         bin = (upper_bound_x - lower_bound_x) / 50  


                        cleanup_time_data <- peak_stat_df %>% dplyr::filter(t_rise > 1 & !is.na(t_rise)) %>% 
                           										dplyr::filter(t_half > 1 & !is.na(t_half)) %>%
                           										dplyr::filter(t_decay > 1 & !is.na(t_decay)) %>%
                           										dplyr::filter(t_half <= 125)
                        time_data_medians<- cleanup_time_data %>%
                        					group_by(Ca, segmentation) %>%
                        					summarise(median_t_rise = median(t_rise,na.rm=TRUE),
											  			median_t_half = median(t_half,na.rm=TRUE),
											  			median_t_decay = median(t_decay,na.rm=TRUE) )
		                           
                        peak_stat_df$segmentation<-factor(peak_stat_df$segmentation, levels=c("marker","activity")) 				
                        cleanup_time_data$segmentation<-factor(cleanup_time_data$segmentation, levels=c("marker","activity")) 				


                        my_medians <- peak_stat_df %>%
									  	group_by(Ca, segmentation) %>%
									  	dplyr::filter(t_half <= 125) %>%
									  summarize(median_amplitude = median(amplitude),
									  			median_tau = median(tau_decay_ms),
									  			median_t_rise = median(t_rise,na.rm=TRUE),
											  			median_t_half = median(t_half,na.rm=TRUE),
											  			median_t_decay = median(t_decay,na.rm=TRUE) ) #,

					   histo_amplitude<-ggplot(peak_stat_df, aes(x=amplitude, group=segmentation, fill=segmentation,colour=segmentation))+
					   					   geom_vline(data=my_medians, aes(xintercept = median_amplitude, colour=segmentation),size = scalefactor * 2, lty = "solid",show_guide = FALSE,alpha=0.9) +
    										   geom_freqpoly(aes(y = ..ncount..), binwidth = bin, alpha = 0.8, boundary = 0, size = scalefactor * 1,show_guide = FALSE) +
    										   geom_histogram(aes(y = ..ncount..), alpha=0.4,binwidth = bin, boundary = 0, position = "identity", colour = "black", size = scalefactor * 0.5,show_guide = TRUE) +
    										   geom_text(aes(label = Ca), x = upper_bound_x+(0.05*upper_bound_x), colour="black",y = 0.6, size=7,hjust=1,check_overlap = TRUE,parse=TRUE,show_guide=FALSE)+
    
					                                        #    geom_vline(data=my_medians, aes(xintercept = median_amplitude, colour=segmentation), size=scalefactor*2, lty="solid" )+
                                        						
					                                         #   geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					#geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression("Peak "*Delta*"F/F"[y]),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "E"
					                                                    )+
					                                            
					                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'), guide= "none")}+
					                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(values = c("activity" = "red", 'marker' = 'black'))}+
					                                            #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.6, 'marker' = 0.4))}+
					                                            
					                                            
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

					     amplitude_histo <<- histo_amplitude




					       lower_bound_x = 0
	                      upper_bound_x = 150
      
                         bin = (upper_bound_x - lower_bound_x) / 50  

					     histo_tau_decay_ms<-ggplot(peak_stat_df, aes(x=tau_decay_ms, group=segmentation, fill=segmentation,colour=segmentation))+
					                                        geom_vline(data=my_medians, aes(xintercept = median_tau, colour=segmentation),size = scalefactor * 2, lty = "solid",show_guide = FALSE,alpha=0.9) +
    										   geom_freqpoly(aes(y = ..ncount..), binwidth = bin, alpha = 0.8, boundary = 0, size = scalefactor * 1,show_guide = FALSE) +
    										   geom_histogram(aes(y = ..ncount..), alpha=0.4,binwidth = bin, boundary = 0, position = "identity", colour = "black", size = scalefactor * 0.5,show_guide = TRUE) +
										#     geom_vline(data=my_medians, aes(xintercept = median_tau, colour=segmentation), size=scalefactor*2, lty="solid" )+
                                        						
					                                        #     geom_histogram(aes(y=..ncount..,alpha=segmentation),binwidth = bin,alpha=0.4, boundary = 0,position="identity",colour="black",size=scalefactor*0.5)+
                                             					# geom_freqpoly(aes(y=..ncount..),binwidth = bin, alpha=0.7,boundary =0,size=scalefactor*1)+
                                             					 
					                                            labs( x=expression(tau[decay]*" (ms)"),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#""expression(tau[decay]*" (ms)"),
					                                                    fill="Segmented by:",
					                                                    tag = "F"
					                                                    )+
					                                            
					                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'), guide= "none")}+
					                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(values = c("activity" = "red", 'marker' = 'black'))}+
					                                            #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.6, 'marker' = 0.4))}+
					                                            
					                                            
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

					     tau_histo <<- histo_tau_decay_ms

			



					       lower_bound_x = 0
	                      upper_bound_x = 150
      
                         bin = (upper_bound_x - lower_bound_x) / 50  



					     histo_t_rise<-ggplot(cleanup_time_data, aes(x=t_rise, group=segmentation, fill=segmentation,colour=segmentation))+
					                                            geom_vline(data=time_data_medians, aes(xintercept = median_t_rise, colour=segmentation), size=scalefactor*2, lty="solid" )+
                                        						
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
					                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'), guide= "none")}+
					                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(values = c("activity" = "red", 'marker' = 'grey21'))}+
					                                            {if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
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

					     histo_t_decay<-ggplot(cleanup_time_data, aes(x=t_decay, group=segmentation, fill=segmentation,colour=segmentation))+
					                                            geom_vline(data=time_data_medians, aes(xintercept = median_t_decay, colour=segmentation), size=scalefactor*2, lty="solid" )+
                                        						
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
					                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'), guide= "none")}+
					                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(values = c("activity" = "red", 'marker' = 'grey21'))}+
					                                            {if(color_switch)scale_alpha_manual(values = c("activity" = 0.5, 'marker' = 0.3))}+
					                                            
					                                            
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
	                      upper_bound_x = 120
      
                         bin = (upper_bound_x - lower_bound_x) / 50 

					     histo_t_half<-ggplot(cleanup_time_data, aes(x=t_half, group=segmentation, fill=segmentation,colour=segmentation))+
					                                        geom_vline(data=my_medians, aes(xintercept = median_t_half, colour=segmentation),size = scalefactor * 2, lty = "solid",show_guide = FALSE,alpha=0.9) +
    										   geom_freqpoly(aes(y = ..ncount..), binwidth = bin, alpha = 0.8, boundary = 0, size = scalefactor * 1,show_guide = FALSE) +
    										   geom_histogram(aes(y = ..ncount..), alpha=0.4,binwidth = bin, boundary = 0, position = "identity", colour = "black", size = scalefactor * 0.5,show_guide = TRUE) +
										#stat_ecdf(size=2,alpha=1, lineend="round")+
					                                             
					                                            labs( x=expression(t[1/2*"y"]),
					                                                    y=expression(N/N[max]),
					                                                    #title="",#"Peak amplitude",
					                                                    fill="Segmented by:",
					                                                    tag = "G"
					                                                    )+
					                                            
					                                          {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'), guide= "none")}+
					                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
					                                            {if(color_switch)scale_fill_manual(values = c("activity" = "red", 'marker' = 'black'))}+
					                                            #{if(color_switch)scale_alpha_manual(values = c("activity" = 0.6, 'marker' = 0.4))}+
					                                            
					                                            
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
					                                            #legend.position = c(.95, .95),
					                                                  #legend.justification = c("right", "bottom"))#,
					                                                   # legend.margin = margin(6, 6, 6, 6))
					                       
					     ggsave(filename=paste0("histo_t_half.jpeg"),plot=histo_t_half, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)

					     thalf_histo <<- histo_t_half

####### SCATTERPLOTS ####### 




                        plot_avgs <- peak_stat_df %>%
									  	group_by(Ca, segmentation) %>%
									  summarize(n = n(),
									         mean_amplitude = mean(amplitude,na.rm=TRUE),
									         sd_amplitude = sd(amplitude,na.rm=TRUE),
									         se_amplitude = sd_amplitude / sqrt(n),
									         
									         mean_tau = mean(tau_decay_ms,na.rm=TRUE),
									         sd_tau = sd(tau_decay_ms,na.rm=TRUE),
									         se_tau = sd_tau / sqrt(n))
					   print(plot_avgs)

						scatter_amplitude_tau<-ggplot(peak_stat_df, aes(x=amplitude, y=tau_decay_ms, group=segmentation, colour=segmentation))+#linetype=segmentation,
						                                            geom_point(size=.5, alpha=0.2)+
						                                            geom_errorbar(data=plot_avgs, aes(x=mean_amplitude,y=mean_tau, ymin = mean_tau - se_tau, ymax = mean_tau + se_tau, group=segmentation, colour=segmentation, width=0.15), alpha=1)+
						                                            geom_errorbarh(data=plot_avgs, aes(x=mean_amplitude, y=mean_tau, xmin = mean_amplitude - se_tau, xmax = mean_amplitude + se_tau, colour=segmentation,height=5), alpha=1)+
						                                            
						                                            geom_point(data=plot_avgs, aes(x=mean_amplitude, y=mean_tau, group=segmentation, colour=segmentation),size=2, alpha=1)+
						                                            labs( y=expression(tau[decay]*" (ms)"),
						                                                    x=expression("Peak "*Delta*"F/F"),
						                                                    title="",#expression(Delta*"t (ms)"),
						                                                    colour="Segmented by:"
						                                                    )+
						                                            
						                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                                            {if(color_switch)scale_colour_manual(values = c("activity" = "red", 'marker' = 'black'))}+
						                                            #scale_linetype_manual(name = "", values = c("activity" = "dashed", 'marker' = 'solid'))+
						                                            coord_cartesian(ylim=c(0,200),xlim=c(0,6))+
						                                            #scale_y_continuous(breaks=c(0,0.5,1.0))+
						                                            #scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
						                            
						                                            theme_tufte()+
						                                            facet_grid(Ca~., labeller = label_parsed)+
						                                            my.theme+
						                                            theme(legend.title = element_text(colour="black", size=28, family="sans"),
						                                                  legend.position = "right",
					                                            	      strip.text.y.right = element_text(size=22, angle = 0))#,
						                                            	  #strip.text.y.right = element_text(angle = 0))
						                                            #legend.position = c(.95, .95),
						                                                  #legend.justification = c("right", "bottom"))#,
						                                                   # legend.margin = margin(6, 6, 6, 6))
						                       
							ggsave(filename=paste0("scatter_amplitude_tau.jpeg"),plot=scatter_amplitude_tau, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_width,height=cumdist_height)








# Goal of this script is to accept a summary data.frame, peakStats, and plot data visualizations for Figure 2 of the optical physiology manuscript. 
library(ggh4x)
library(ggpubr)
library(ggforce)


##### AVG PPR WAVEFORMS ######
						clean_peaks_normTime<- clean_peaks_normTime %>% dplyr::filter(new_normTime >= -0.055, new_normTime <= 0.30) 
						clean_peaks_normTime$segmentation <- factor(clean_peaks_normTime$segmentation, levels = rev(segmentation_levels))

						
						avgtrace_Plot<-ggplot(clean_peaks_normTime, aes(x=new_normTime, y = dFF,  group=segmentation, colour=segmentation, linetype = segmentation))+
						                            #geom_path(colour="grey61",size=0.4,alpha=0.2)+
						                            stat_summary_bin(aes(group=segmentation), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.8)+
						                            #geom_vline(xintercept = 0, lty="dashed",colour='grey21',size=0.9)+
						                            
						                           	labs( x="Normalized time (s)",
						                                  y=expression(Delta*"F/F"),
						                                  #colour="Segmented by:",
					                                                    tag = "D")+
						                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						                            {if(color_switch)scale_colour_manual(values=c("activity" = "red", 'marker' = 'black'))}+
						                            scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						                            scale_x_continuous(breaks = c(0,0.3),labels=c("0","0.3"),limits=c(-0.15,0.6))+
						                            theme_tufte()+
						                            #guides(colour="none",fill="none")+
						                            coord_cartesian(ylim=c(-0.05,0.55),xlim=c(-0.1,0.4))+
						                            my.theme+
						                            facet_nested_wrap(~Ca + segmentation, nrow = 1, ncol=8, labeller = labeller(Ca = label_parsed, 
						                            									        segmentation = function(x) {rep("", length(x))}))+
						                            theme(legend.title = element_text(colour="black", size=28, family="sans"),
						                            	legend.position="none")

						                            #avg_peaks  <<- avg_PP_tracePlot
						                             #       rm(avg_PP_tracePlot
						                                

						                             ggsave(filename=paste0("averagePeaks_.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)
                         

						avg_traces <<- avgtrace_Plot
					     

						# clean_peaks_long<- clean_peaks_normTime %>% dplyr::filter(tau_decay_ms >200)




						 # avgtrace_Plot<-ggplot(clean_peaks_long, aes(x=new_normTime, y = dFF,  group=segmentation, colour=segmentation, linetype = segmentation))+
						 #                            #geom_path(colour="grey61",size=0.4,alpha=0.2)+
						 #                            stat_summary_bin(aes(group=segmentation), geom="line", fun=mean, binwidth=interFrame.var,size=2,alpha=0.8)+
						 #                            geom_vline(xintercept = 0, lty="dashed",colour='black',size=0.7)+
						                            
						 #                           	labs( x="Normalized time (s)",
						 #                                  y=expression(Delta*"F/F"),
						 #                                  colour="Segmented by:")+
						 #                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						 #                            {if(color_switch)scale_colour_manual(values=c("activity" = "red", 'marker' = 'black'))}+
						 #                            scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						 #                            scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
						 #                            theme_tufte()+
						 #                            #guides(colour="none",fill="none")+
						 #                            #coord_cartesian(ylim=c(-0.1,2.0),xlim=c(-1.1,1.1))+
						 #                            my.theme+
						 #                            facet_grid(~Ca, labeller = label_parsed)+
						 #                            theme(legend.title = element_text(colour="black", size=28, family="sans"))

						 #                            #avg_peaks  <<- avg_PP_tracePlot
						 #                             #       rm(avg_PP_tracePlot)
						                                

						 #                             ggsave(filename=paste0("averagePeaks_long.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)

						   



						 #    clean_peaks_short<- clean_peaks_normTime %>% dplyr::filter(tau_decay_ms <30)





						 #   avgtrace_Plot<-ggplot(clean_peaks_short, aes(x=new_normTime, y = dFF,  group=segmentation, colour=segmentation, linetype = segmentation))+
						 #                            #geom_path(colour="grey61",size=0.4,alpha=0.2)+
						 #                            stat_summary_bin(aes(group=segmentation), geom="line", fun=mean, binwidth=interFrame.var,size=2,alpha=0.8)+
						 #                            geom_vline(xintercept = 0, lty="dashed",colour='black',size=0.7)+
						                            
						 #                           	labs( x="Normalized time (s)",
						 #                                  y=expression(Delta*"F/F"),
						 #                                  colour="Segmented by:")+
						 #                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
						 #                            {if(color_switch)scale_colour_manual(values=c("activity" = "red", 'marker' = 'black'))}+
						 #                            scale_linetype_manual(name = "", values = c("activity" = "solid", 'marker' = 'solid'),guide="none")+
						 #                            scale_x_continuous(breaks = c(0,0.5),limits=c(-0.15,0.6))+
						 #                            theme_tufte()+
						 #                            #guides(colour="none",fill="none")+
						 #                            #coord_cartesian(ylim=c(-0.1,2.0),xlim=c(-1.1,1.1))+
						 #                            my.theme+
						 #                            facet_grid(~Ca, labeller = label_parsed)+
						 #                            theme(legend.title = element_text(colour="black", size=28, family="sans"))

						 #                            #avg_peaks  <<- avg_PP_tracePlot
						 #                             #       rm(avg_PP_tracePlot)
						                                

						 #                             ggsave(filename=paste0("averagePeaks_short.jpeg"),plot=avgtrace_Plot, device="jpeg",dpi=600, bg="white",units="in",width=cumdist_height,height=cumdist_width)

list.output = list(ROI_totals = getROIs, ROI_probability = ROIs_with_peaks, ROIs_per_vid_summary = check_ROIs_per_vid, Medians = my_medians, time_Medians = time_data_medians, Averages = plot_avgs, Original_stats =peak_stat_df)
list.output
}


         