subDir <- "single_ROI_display_v9"   


source(paste0(path,"sub-functions/setFigureDirectory.R")) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"peakAnalysis/png_plotter_v5_ROIoverlay_syntransmission.R"))
#source(paste0(path,"peakAnalysis/plot_dualRecords_fxn.R"))
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))



library(ggforce)
library(scales)
library(ggthemes)
library(ggpmisc)
library(ggpubr)



single_ROI_display<- function(ROI_traces, interSpike_thresh = NULL, peak_groupers,groupers, positive_ROI_groupers,plotBy, levels, secondAxis = NULL, color_override=color_override,comparisons, save_ROIs = ROIs_to_save,tag_var = tmp_tag, display_transmission = disp_trans){

									window_boundary = c(0.01,0.25)
									sd_threshold = 3.5

									minpositive = function(x) min(x[x > 0 ])
									cumplus <- function(y) Reduce(function(a,b) a + b > 0, y, 0, accum=TRUE)[-1]
			

										 vid_keys_to_save = unique(str_extract(save_ROIs, "dish(.*?)region\\d\\d") )
				                    print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))
				                    if (!is.null(interSpike_thresh)) { interSpike_thresh = interSpike_thresh
				                                                        } else {
				                                                            interSpike_thresh = 0.25
				                                                            } ## set threshold for accepting interSpike interval as transmission
				                    secondAxis_switch = !is.null(secondAxis)
				                    color_switch = !is.null(color_override)
				                    comparison_switch = !is.null(comparisons)
				                    

				                    ## generate the keys for tracking ROIs
				                    #barplot_width = 16
				                    #lineplot_width = 16
				                    scalefactor = 0.75
				                    plot_dim = 16

				                    total_ROIs = length(save_ROIs)
				                    tag_var_list = c("B")

				                    for(i in 1:total_ROIs) {
					                            current_tag = tag_var_list[i]
					                            subset_traces = ROI_traces %>% dplyr::filter(trackROI_key %in% save_ROIs[i])
					                            if(unique(subset_traces$vid_key) == 'dish7-plate07-region04'){

					                            	traces_df_fix_brokenROI<- subset_traces %>% dplyr::filter(chemical_condition == "APV", vid_key == "dish7-plate07-region04", absoluteTime >= 6) ###this dataset had a focus issue that is messing with quantification
									                traces_df_revise_brokenROI <- subset_traces %>% dplyr::filter(chemical_condition == "ctrl",vid_key == "dish7-plate07-region04")
									                #traces_df_omit_brokenROI<- traces_df %>% dplyr::filter(vid_key != "dish7-plate07-region04" )
													traces_df_fixed <- bind_rows(traces_df_fix_brokenROI,traces_df_revise_brokenROI)

	                            					subset_traces <- traces_df_fixed %>% mutate(new_absoluteTime = absoluteTime - min(absoluteTime)) %>% select(-absoluteTime)
	                            					names(subset_traces)[names(subset_traces) == "new_absoluteTime"] <- "absoluteTime" 
	                            					rm(traces_df_fixed,traces_df_revise_brokenROI,traces_df_fix_brokenROI)
					                            }

					                            # if(i == 1){
					                            # 	facet_switch = TRUE
					                            # } else {
					                            # 	facet_switch = FALSE
					                            # }

					                            	#intensity_notPeak = ifelse(peak_idx == "notPeak",intensity_smooth, NA),
      												#baseline_v0 = na.locf(intensity_notPeak, fromLast = FALSE, na.rm = FALSE),
      												#0baseline = rollapply(baseline_v0, window, mean, fill = NA, align="center")


      												window = 0.75/mean(subset_traces$interFrame,na.rm=TRUE)

															minpositive = function(x) min(x[x > 0 ])



															        tmp_trace <- subset_traces %>% group_by(sensor,chemical_condition,ROINumber)
															        print(tmp_trace)

															        get_baseline_indices <- tmp_trace %>% dplyr::filter(sensor == "JF646",alpha_str == "notPeak")
															        baseline_times <- get_baseline_indices %>% ungroup() %>% select(chemical_condition, trackROI_key, absoluteTime)
															        baseline_pairs<-inner_join(tmp_trace,baseline_times)
															        #baseline_pairs<- tmp_trace %>% dplyr::filter(absoluteTime %in% baseline_times)
															        get_hline<- baseline_pairs %>% dplyr::filter(sensor == "JF646", chemical_condition == "ctrl") %>% summarise(hline_mean = mean(dFF,na.rm=TRUE),
															        																											hline_sd = sd_threshold*sd(dFF,na.rm=TRUE),
															        																											hline_tresh = hline_mean+hline_sd)
															        hline_val <- get_hline$hline_tresh
															        hline_mean<- get_hline$hline_mean


															        print("Checking our hline approach.")
															        #print(baseline_times)
															        #print(unique(baseline_pairs$sensor))
															        #print(get_hline)
															        print(hline_val)

																	#baseline_intensity <- round(mean(extract_baseline_intensity$dFF,na.rm=TRUE),2)  															        
															        #print(tmp_trace)
															        								
															        graphTitle = first(unique(tmp_trace$trackROI_key) )
															        print(paste0("Current ROI is: ", graphTitle))


															                    def_Glu_window_boundaries <- tmp_trace %>% 
															                                      dplyr::filter(peakID != "NotPeak") %>%
															                                      dplyr::filter(sensor == "GluSnFR3") %>%
															                                      group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
															                                      slice(which.max(dFF)) %>% 
															                                      ungroup() %>% 
															                                      group_by(chemical_condition,trackROI_key) %>% 
															                                      arrange(absoluteTime) %>%
															                                      mutate(#Glu_peakID = paste0("Glu_",peakID),
															                                      		 check_time_start = absoluteTime-window_boundary[1],
															                                      		 check_time_end = absoluteTime + window_boundary[2]#interSpike_thresh
															                                               ) %>%
															                                      ungroup() %>%
															                                      select(sensor, chemical_condition,trackROI_key,absoluteTime,dFF,check_time_start,check_time_end) #Glu_peakID,

															                    Ca_trace <- tmp_trace %>% dplyr::filter(sensor == "JF646")

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
															                    selected_intervals <- list()
															                    ## use for loop (bleh) to define bins of data.frame
															                    for (i in 1:n_intervals) {

															                    		current_chem_condition = start_times$chemical_condition[i]
															                    		bin_number = i
															                    		bin_start = start_times$check_time_start[i]
															                    		bin_end = end_times$check_time_end[i]
															                    		selected_data <- Ca_trace %>% ungroup() %>% 
															                    										dplyr::filter(chemical_condition == current_chem_condition, absoluteTime >= bin_start & absoluteTime <= bin_end) %>%
															                    										mutate(window_start = first(absoluteTime),
															                    												window_end = last(absoluteTime),
															                    												ribbon_window = 1,
															                    												bin = bin_number)

															                    		selected_intervals[[i]] <- selected_data

															                    }

															                   	Ca_trace_binned <- bind_rows(selected_intervals) 

															                   	Ca_plot_elements<- Ca_trace_binned %>% ungroup() %>% select(sensor, chemical_condition, trackROI_key,window_start, window_end)
															                   	Ca_plot_elements<-unique(Ca_plot_elements)

															                   	get_JF_peakPositions <- Ca_trace_binned %>% 
															                   						group_by(sensor,chemical_condition,trackROI_key,bin) %>%
															                   						slice(which.max(dFF)) %>%
															                   						#dplyr::filter(dFF >= 0.05) %>%
															                   						mutate(peak_color = ifelse(dFF >= hline_val, "blue","white"))#,
															                   								#peak_shape = ifelse(dFF >= 0.1, 21, 4))





															                   	if(graphTitle == "dish7-plate07-region04-ROI0019") {
															                   					graphSubtitle_trace = ""#"Single Spine Glutamate / Calcium Dynamics"
															                   					graphSubtitle_amp = "Single Spine Dose-Response"
															                   	} else if(graphTitle == "dish8-plate07-region05-ROI0020") {
															                   		graphSubtitle_trace = "Control Wash-in"
															                   		graphSubtitle_amp = ""
															                   	} else if(graphTitle == "dish7-plate10-region04-ROI0016") {
															                   		graphSubtitle_trace = "AP5 Wash-in"
															                   		graphSubtitle_amp = ""
															                   	} else {
																					graphSubtitle_trace = ""
															                   		graphSubtitle_amp = ""

															                   	}
															                	

															                	#define baseline dataframe

															                   	Glu_baseline<- baseline_pairs %>% dplyr::filter(sensor == "GluSnFR3") %>% ungroup() %>% select(chemical_condition, trackROI_key,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
															                    colnames(Glu_baseline)[colnames(Glu_baseline) == "dFF"] <- "Glu_dFF"

															                    JF_baseline<- baseline_pairs %>% dplyr::filter(sensor == "JF646") %>% ungroup() %>% select(chemical_condition, trackROI_key,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
															                    colnames(JF_baseline)[colnames(JF_baseline) == "dFF"] <- "JF_dFF"

															                    check_baseline_dist<- left_join(Glu_baseline,JF_baseline) %>% mutate(peak_color = "baseline")


															                    ## define peak position dataframe
															                    get_Glu_peakPositions = def_Glu_window_boundaries              


															                    Glu_peaks <- get_Glu_peakPositions %>% ungroup() %>% select(chemical_condition, trackROI_key,dFF) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
															                    colnames(Glu_peaks)[colnames(Glu_peaks) == "dFF"] <- "Glu_dFF"
															                    JF_peaks <- get_JF_peakPositions %>% ungroup() %>% select(chemical_condition, trackROI_key,dFF,peak_color) %>% group_by(chemical_condition, trackROI_key) %>% mutate(index = row_number())
															                    colnames(JF_peaks)[colnames(JF_peaks) == "dFF"] <- "JF_dFF"

															                    check_amplitude_dist = left_join(Glu_peaks,JF_peaks, by=c("index","chemical_condition","trackROI_key"))
															                    check_amplitude_dist = bind_rows(check_amplitude_dist,check_baseline_dist)
															                    check_amplitude_dist$chemical_condition = factor(check_amplitude_dist$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															
															                     if(disp_trans == FALSE) {
																					 	strip_text_arg_y = element_text(colour="black", size=34*scalefactor, family="sans", angle=0)
																					 	strip_text_arg_x = element_blank()
																					 	legend_arg = "right"
																					 	break_arg_x = c(0,0.5,1,1.5)
																					 	break_arg_y = c(0,0.4,0.8)
																					 	xlim_val = 1.8
																					 	ylim_val = 1

																					 } else if (disp_trans == TRUE) {
																					 	strip_text_arg_y = element_blank()
																					 	strip_text_arg_x = element_text(colour="black", size=34*scalefactor, family="sans")
																					 	legend_arg = "none"
																					 	break_arg_x = c(0,0.5,1,1.5)
																					 	xlim_val = 1.7
																					 	ylim_val = 0.8
																					 	break_arg_y = c(0,0.4,0.8)
																					 	

																					 }

																				
															                    amplitude_scatterplot<-ggplot(check_amplitude_dist, aes(x=Glu_dFF, y=JF_dFF, group=interaction(trackROI_key,chemical_condition),colour=peak_color,alpha=peak_color))+
																	                                                geom_smooth(data=subset(check_amplitude_dist,peak_color == "blue"), method = "lm",formula = y~ x, se = FALSE, colour = "black", size=1,alpha=0.3)+
																	                                                geom_point(size=4)+
																	                                                {if(disp_trans == FALSE)stat_cor( color = "black", geom = "text",label.x = 0.15*xlim_val,label.y=0.95,size=7, p.accuracy = 0.001)}+ #aes(label = after_stat(rr.label)),
																	                                                #{if(disp_trans == FALSE)stat_regline_equation(label.x = 0.15*xlim_val, label.y = 0.85,color="black",size=7)} +
																	                                                geom_hline(yintercept =hline_val,lty="dashed",size=0.5)+
																	                                                {if(disp_trans == FALSE)geom_label(label = expression("3.5" %*% sigma[baseline]), 
																								                               y = hline_val,	 x = 0.85*xlim_val, size = 7, colour = 'black',fontface='italic', fill="white", label.size=NA)}+
																	                                      labs( x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
																	                                                  y=expression(JF[646]*" amplitude ("*Delta*"F/F)"),
																	                                                  ##title="",
																	                                                  subtitle=graphSubtitle_amp,#"Single Spine Response",#paste0("Spine",gsub("ROI00", "  #", graphTitle)),
																	                                                  tag=tag_var[2])+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
																	                                           scale_colour_manual(name = "", breaks=c("blue","white", "baseline"),
																	                                           									values = c("blue"='blue',"white" = "grey61","baseline" = "grey61"),
																	                                           									labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"))+
																	                                           scale_alpha_manual(name = "", values = c("blue"=0.9,"white" = 0.4,"baseline" = 0.4),labels=c("blue" = "Transmission", "white" = "No Transmission","baseline" = "Baseline"), guide="none")+																	                                           
																	                                           coord_cartesian(xlim=c(-0.1,xlim_val),ylim=c(-0.1,ylim_val))+
																	                                           scale_y_continuous(breaks=break_arg_y)+#c(0,0.2,0.4,0.6,0.8,1.0))+
																	                                           scale_x_continuous(breaks=break_arg_x)+#scale_y_continuous(limits=c())
																	                                           facet_grid(~chemical_condition,labeller = label_parsed)+
																	                                           theme_tufte()+
																	                                           my.theme +
																	                                           theme(legend.position = legend_arg,
																	                                           		 strip.text.x = strip_text_arg_x)
																	                                # save_plot(filename=paste0("amplitude_corr_",graphTitle),plot=amplitude_scatterplot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)



															            

															min_val = -0.1
															
															facet.labels <- c("iGluSnFR3","JF646-BAPTA-AM")
															names(facet.labels)<- c("GluSnFR3", "JF646")
															tmp_trace$chemical_condition = factor(tmp_trace$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															Ca_plot_elements$chemical_condition = factor(Ca_plot_elements$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															get_Glu_peakPositions$chemical_condition = factor(get_Glu_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															get_JF_peakPositions$chemical_condition = factor(get_JF_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", expression("100"~mu*"M AP5")))
															hline_df<- data.frame(sensor = "JF646", hline = hline_val, min = min_val)

															get_times<- tmp_trace %>% dplyr::filter(sensor == 'JF646') %>% ungroup() %>% select(sensor,chemical_condition,trackROI_key,absoluteTime) %>% 
																															group_by(sensor, chemical_condition,trackROI_key) %>%
																															summarise(max_time = max(absoluteTime,na.rm=TRUE))
															hline_df<- left_join(get_times,hline_df)





															 custom_y <- tribble(
															  ~sensor, ~absoluteTime, ~dFF,
															  "GluSnFR3", 0, min_val,
															  "GluSnFR3", 0, 1.2,
															  "JF646", 0, min_val,
															  "JF646", 0, 0.7
															                    )

															 my_y_breaks <- function(x) {if(max(x) <= 1.1)  seq(0.0,1.0,by=0.3) else if(max(x) <= 2 & max(x) >= 1.2) seq(0.0,2,by=0.5) else seq(0.0,7.5,by=2.5)}
															 my_breaks <- function(x) { if(max(x) < 20)  seq(0,15,5) else seq(0,25,5) }
															 formatter <- function(...){
																						  function(x) format(round(x, 1), ...)
																						}


															 tmp_trace$sensor = as.factor(as.character(tmp_trace$sensor))
															 levels(tmp_trace$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA"))
															 
															Ca_plot_elements$sensor = as.factor(as.character(Ca_plot_elements$sensor))
															 levels(Ca_plot_elements$sensor) <- c(expression(JF[646]*"-BAPTA"))

															 get_Glu_peakPositions$sensor = as.factor(as.character(get_Glu_peakPositions$sensor))
															 levels(get_Glu_peakPositions$sensor) <- c("iGluSnFR3")

															 get_JF_peakPositions$sensor = as.factor(as.character(get_JF_peakPositions$sensor))
															 levels(get_JF_peakPositions$sensor) <- c(expression(JF[646]*"-BAPTA"))


															 custom_y$sensor = as.factor(as.character(custom_y$sensor))
															 levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA"))

															 hline_df$sensor = as.factor(as.character(hline_df$sensor))
															 levels(hline_df$sensor) <- c(expression(JF[646]*"-BAPTA"))
															 

															 if(disp_trans == FALSE) {
															 	strip_text_arg_y = element_blank()#element_text(colour="black", size=34*scalefactor, family="sans", angle=0)
															 	strip_text_arg_x = element_blank()
															 	legend_arg = "right"
															 	legend_just = "left"
															 } else if (disp_trans == TRUE) {
															 	strip_text_arg_y = element_blank()
															 	strip_text_arg_x = element_text(colour="black", size=34*scalefactor, family="sans")
															 	legend_arg = "none"
															 	legend_just = "left"
															 }

															tracePlot<-ggplot(tmp_trace, aes(x=absoluteTime, y = dFF, group = sensor)) +
															            #geom_rect(data=hline_df,aes(xmin=-Inf, xmax=Inf,ymin = -Inf, ymax=hline_df$hline), fill="grey71",alpha=0.2,inherit.aes=FALSE)+
															            #geom_hline(yintercept = hline_mean, lty="solid", colour='black',size=1.1)+
															            geom_rect(data=Ca_plot_elements, aes(xmin = window_start, xmax = window_end,group=chemical_condition), ymin = -Inf, ymax = Inf,alpha=0.14, fill="grey21",inherit.aes=FALSE)+
															            geom_vline(data=Ca_plot_elements,aes(xintercept = window_start,group=chemical_condition), colour = "grey21", alpha = 0.14, linetype = "solid", size =0.1) +
															            geom_vline(data=Ca_plot_elements,aes(xintercept = window_end,group=chemical_condition), colour = "grey21", alpha = 0.14, linetype = "solid", size =0.1) + 
															            #geom_vline(data=sliceTime, aes(xintercept = absoluteTime), colour = "grey21", alpha = 0.5, linetype = "dashed", size =0.8) +
															            #{if(display_transmission)geom_text(data=sliceGluPeaks_v2,aes(label = "X",x = absoluteTime, y = yposX), colour="red", size= 8,alpha=0.75)}+
															            #{if(display_transmission)geom_text(data=sliceGluPeaks_v2,aes(label = "O",x = absoluteTime, y = yposO), colour="black", size= 8,alpha=0.75)}+
															            geom_blank(data=custom_y,aes(absoluteTime,dFF))+
															            geom_line(aes(colour = sensor),size = 1.2) + 
															           
															            geom_point(data=get_Glu_peakPositions, aes(x=absoluteTime,y=dFF), shape=23, fill="blue",alpha=0.7,colour="black", size=4)+
															            geom_point(data=get_JF_peakPositions, aes(x=absoluteTime,y=dFF,fill=peak_color,shape=peak_color), colour="black",alpha=0.7, size=4)+ #colour="blue",fill="blue", shape=23,
															            
															            labs( x="Time (s)",
															                  y=expression(Delta*"F/F"), 
															                  tag = tag_var[1],
															                  subtitle = graphSubtitle_trace)+#,
															                  #subtitle=paste0("Spine",gsub("ROI00", "  #", graphTitle)))+
															            scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA")))+
															            scale_fill_manual(values = c("blue","grey71"),labels=c("blue","white"),guide='none')+
															            scale_shape_manual(values = c(23,23),labels=c("blue","white"), guide='none')+
															            scale_x_continuous(expand=c(0,0))+
															            coord_cartesian(clip='off')+
															            theme_tufte()+
															            my.theme+
															            
															             theme(legend.position=legend_arg,
															             		strip.text.x=strip_text_arg_x,
															             		strip.text.y= strip_text_arg_y,
															             		legend.justification = legend_just)+
															            facet_grid(sensor~chemical_condition, labeller = label_parsed,scales="free_y")+
															            scale_y_continuous(labels = formatter(nsmall=1), breaks= my_y_breaks) #expand=c(0,0),
															            #+
															             


															#save_plot(filename=paste0(graphTitle,"dualRecording_windowed_maxima_back"),plot=tracePlot, plot_width=plot_dim*2.25,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)
															


															gs = list(tracePlot,  #1
										                                 amplitude_scatterplot)  #3
										                        margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))              
										                    



										                        hlay <- rbind(c(1,1,1,1,1,1,1,2,2,2,2),
										                                       c(1,1,1,1,1,1,1,2,2,2,2),
										                                       c(1,1,1,1,1,1,1,2,2,2,2),
										                                       c(1,1,1,1,1,1,1,2,2,2,2)
										                                   )
										                
										                    

										                    
										                        ROI_report = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
										                        save_plot(filename=paste0("ROI_report_",graphTitle,"_"),plot=ROI_report, plot_width=plot_dim*3,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 200)
													



					}
output<-list(tracePlot,amplitude_scatterplot)
output
#output <- list(Glu_peaks,JF_peaks,check_amplitude_dist)
#output 
#tmp_trace
#Ca_plot_elements#window_boundaries
}



# interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
# 															                                             isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
# 															                                             timeCheck = ifelse(isValid == TRUE & interSpike < interSpike_thresh, TRUE, FALSE),
# 															                                              isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
# 															                                                                          sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
# 															                                                                          sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
# 															                                              checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)
# 															                                             




															                    # checkInterspike<- slicePeaks %>% dplyr::filter(timeCheck == TRUE)
															                    # getOrphans<- slicePeaks %>% dplyr::filter(timeCheck == FALSE)

															                    # sliceJF646_peaks<- slicePeaks %>% dplyr::filter(sensor == "JF646")
															                    # get_JF_peakPositions = sliceJF646_peaks
    # sliceGluPeaks_v2<- tmp_trace %>%  dplyr::filter(peakID != "NotPeak") %>% 
	# 														                                                          group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
	# 														                                                          slice(which.max(dFF)) %>% 
	# 														                                                          ungroup() %>% 
	# 														                                                          group_by(chemical_condition,trackROI_key) %>% 
	# 														                                                          arrange(absoluteTime) %>%
	# 														                                                          mutate(interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
	# 														                                                                    isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
	# 														                                                                     timeCheck = ifelse(isValid == TRUE & interSpike < interSpike_thresh, TRUE, FALSE),
	# 														                                                                      isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
	# 														                                                                                                  sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
	# 														                                                                                                  sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
	# 														                                                                      checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)   ) %>%
	# 														                                                        dplyr::filter(sensor == "GluSnFR3")


	# sliceGluPeaks_v2$timeCheck[is.na(sliceGluPeaks_v2$timeCheck)] <- FALSE
	# 														                    sliceGluPeaks_v2$sensor[sliceGluPeaks_v2$sensor == "GluSnFR3"] <- "JF646"
	# 														                    sliceGluPeaks_v2 <- sliceGluPeaks_v2 %>% ungroup() %>% group_by(chemical_condition) %>%mutate(yposX = 1.0*ifelse(isTransmission == FALSE, 0.6, NA),
	# 														                                                                                yposO = 1.0*ifelse(isTransmission == TRUE, 0.6, NA))

	#sliceTime$chemical_condition = factor(sliceTime$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
															#sliceGluPeaks_v2$chemical_condition = factor(sliceGluPeaks_v2$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
																	#                             if(i == 1){
					# #### save to global env
					#                              traceB  <<- tracePlot

					#                             } else if(i == 2){
					#                              trace_APV_F <<- tracePlot
					                             	
                    # 							} else if(i == 3){

                    # 							}
				 # sliceTime$sensor = as.factor(as.character(sliceTime$sensor))
															 # levels(sliceTime$sensor) <- c(expression(JF[646]*"-BAPTA"))

															 # sliceGluPeaks_v2$sensor = as.factor(as.character(sliceGluPeaks_v2$sensor))
															 # levels(sliceGluPeaks_v2$sensor) <- c(expression(JF[646]*"-BAPTA"))
															 # tmp_baseline_adj$sensor = as.factor(as.character(tmp_baseline_adj$sensor))
															 # levels(tmp_baseline_adj$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA"))