subDir <- "single_ROI_display_v6"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_v5_ROIoverlay_syntransmission.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/plot_dualRecords_fxn.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/save_plot_as_jpeg_and_emf.R")





single_ROI_display<- function(ROI_traces, interSpike_thresh = NULL, peak_groupers,groupers, positive_ROI_groupers,plotBy, levels, secondAxis = NULL, color_override=color_override,comparisons, save_ROIs = ROIs_to_save,tag_var = tmp_tag, display_transmission = disp_trans){


										 vid_keys_to_save = unique(str_extract(save_ROIs, "dish(.*?)region\\d\\d") )
				                    print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))
				                    if (!is.null(interSpike_thresh)) { interSpike_thresh = interSpike_thresh
				                                                        } else {
				                                                            interSpike_thresh = 0.4
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

															        tmp_trace <- subset_traces %>% group_by(sensor,chemical_condition,ROINumber) #%>% 
															        								
															        graphTitle = first(unique(tmp_trace$trackROI_key) )



															                    slicePeaks <- tmp_trace %>% 
															                                      dplyr::filter(peakID != "NotPeak") %>% 
															                                      group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
															                                      slice(which.max(dFF)) %>% 
															                                      ungroup() %>% 
															                                      group_by(chemical_condition,trackROI_key) %>% 
															                                      arrange(absoluteTime) %>%
															                                      mutate(
															                                             interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
															                                             isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
															                                             timeCheck = ifelse(isValid == TRUE & interSpike < interSpike_thresh, TRUE, FALSE),
															                                              isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
															                                                                          sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
															                                                                          sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
															                                              checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)
															                                               )

															                    checkInterspike<- slicePeaks %>% dplyr::filter(timeCheck == TRUE)
															                    getOrphans<- slicePeaks %>% dplyr::filter(timeCheck == FALSE)

															                    sliceJF646_peaks<- slicePeaks %>% dplyr::filter(sensor == "JF646")
															                    get_JF_peakPositions = sliceJF646_peaks


															                    sliceGluPeaks_v2<- tmp_trace %>%  dplyr::filter(peakID != "NotPeak") %>% 
															                                                          group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
															                                                          slice(which.max(dFF)) %>% 
															                                                          ungroup() %>% 
															                                                          group_by(chemical_condition,trackROI_key) %>% 
															                                                          arrange(absoluteTime) %>%
															                                                          mutate(interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
															                                                                    isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
															                                                                     timeCheck = ifelse(isValid == TRUE & interSpike < interSpike_thresh, TRUE, FALSE),
															                                                                      isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
															                                                                                                  sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
															                                                                                                  sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
															                                                                      checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)   ) %>%
															                                                        dplyr::filter(sensor == "GluSnFR3")

															                    get_Glu_peakPositions = sliceGluPeaks_v2                                      
															                    sliceGluPeaks_v2$timeCheck[is.na(sliceGluPeaks_v2$timeCheck)] <- FALSE
															                    sliceGluPeaks_v2$sensor[sliceGluPeaks_v2$sensor == "GluSnFR3"] <- "JF646"
															                    sliceGluPeaks_v2 <- sliceGluPeaks_v2 %>% ungroup() %>% group_by(chemical_condition) %>%mutate(yposX = 1.0*ifelse(isTransmission == FALSE, 0.6, NA),
															                                                                                yposO = 1.0*ifelse(isTransmission == TRUE, 0.6, NA))
															                    sliceTime <- sliceGluPeaks_v2 %>% select(absoluteTime) %>% mutate(sensor = "JF646")


															                    num_transmission_success = sum(sliceGluPeaks_v2$isTransmission)
															                    num_transmission_attempt = nrow(sliceGluPeaks_v2)
															                    transmission_rate = round(num_transmission_success/num_transmission_attempt,2)
															                    manual_transmission_rate = round( (num_transmission_success+7)/num_transmission_attempt,2)



															facet.labels <- c("iGluSnFR3","JF646-BAPTA-AM")
															names(facet.labels)<- c("GluSnFR3", "JF646")
															tmp_trace$chemical_condition = factor(tmp_trace$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
															sliceTime$chemical_condition = factor(sliceTime$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
															sliceGluPeaks_v2$chemical_condition = factor(sliceGluPeaks_v2$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
															get_Glu_peakPositions$chemical_condition = factor(get_Glu_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))
															get_JF_peakPositions$chemical_condition = factor(get_JF_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Control", "AP5"))


															 
															 custom_y <- tribble(
															  ~sensor, ~absoluteTime, ~dFF,
															  "GluSnFR3", 0, -0.1,
															  "GluSnFR3", 0, 1,
															  "JF646", 0, -0.1,
															  "JF646", 0, 0.75
															                    )

															 my_y_breaks <- function(x) {if(max(x) <= 1.1)  seq(0.0,1.0,by=0.5) else if(max(x) <= 5 & max(x) >= 1.2) seq(0.0,5,by=2.5) else seq(0.0,7.5,by=2.5)}
															 my_breaks <- function(x) { if(max(x) < 20)  seq(0,15,5) else seq(0,25,5) }

															 tmp_trace$sensor = as.factor(as.character(tmp_trace$sensor))
															 levels(tmp_trace$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-AM"))
															 
															 sliceTime$sensor = as.factor(as.character(sliceTime$sensor))
															 levels(sliceTime$sensor) <- c(expression(JF[646]*"-BAPTA-AM"))

															 sliceGluPeaks_v2$sensor = as.factor(as.character(sliceGluPeaks_v2$sensor))
															 levels(sliceGluPeaks_v2$sensor) <- c(expression(JF[646]*"-BAPTA-AM"))

															 get_Glu_peakPositions$sensor = as.factor(as.character(get_Glu_peakPositions$sensor))
															 levels(get_Glu_peakPositions$sensor) <- c("iGluSnFR3")

															 get_JF_peakPositions$sensor = as.factor(as.character(get_JF_peakPositions$sensor))
															 levels(get_JF_peakPositions$sensor) <- c(expression(JF[646]*"-BAPTA-AM"))


															 custom_y$sensor = as.factor(as.character(custom_y$sensor))
															 levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-AM"))

															 # tmp_baseline_adj$sensor = as.factor(as.character(tmp_baseline_adj$sensor))
															 # levels(tmp_baseline_adj$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-AM"))
															 if(disp_trans == FALSE) {
															 	strip_text_arg = element_blank()
															 	legend_arg = "right"
															 } else if (disp_trans == TRUE) {
															 	strip_text_arg = element_text(colour="black", size=24, family="sans")
															 	legend_arg = "none"
															 }

															tracePlot<-ggplot(tmp_trace, aes(x=absoluteTime, y = dFF, group = sensor)) +
															            geom_vline(data=sliceTime, aes(xintercept = absoluteTime), colour = "grey21", alpha = 0.5, linetype = "dashed", size =0.8) +
															            {if(display_transmission)geom_text(data=sliceGluPeaks_v2,aes(label = "X",x = absoluteTime, y = yposX), colour="red", size= 8,alpha=0.75)}+
															            {if(display_transmission)geom_text(data=sliceGluPeaks_v2,aes(label = "O",x = absoluteTime, y = yposO), colour="black", size= 8,alpha=0.75)}+
															            geom_blank(data=custom_y,aes(absoluteTime,dFF))+
															            geom_line(aes(colour = sensor),size = 1.2) + 
															           
															            geom_point(data=get_Glu_peakPositions, aes(x=absoluteTime,y=dFF), shape=23, fill="blue",alpha=0.6,colour="blue", size=5)+
															            geom_point(data=get_JF_peakPositions, aes(x=absoluteTime,y=dFF), shape=23, fill="blue",alpha=0.6,colour="blue", size=5)+
															            
															            labs( x="Time (s)",
															                  y=expression(Delta*"F/F"), #"Intensity (a.u.)",#
															                  tag = tag_var#,#current_tag#,
															                  #subtitle=bquote("Estimated P"[transmission]~"="~.(transmission_rate) ~"or"~.(manual_transmission_rate))
															                  )+
															            scale_colour_manual(name = "", values = c( "forestgreen", "#FF33FF"),labels=c(expression("iGluSnFR3"),expression(JF[646]*"-BAPTA-AM")))+
															            theme_tufte()+
															            my.theme+
															            #{if(facet_switch)facet_grid(sensor~., labeller = label_parsed,scales="free_y")}+
															            #facet_grid(sensor~chemical_condition, labeller = label_parsed,scales="free_y")

															             theme(legend.position=legend_arg,
															             		legend.justification="center",
															             		strip.text.x=strip_text_arg,
															             		strip.text.y= element_blank())+#element_text(colour="black", size=30, family="sans",angle=0) )+
															            facet_grid(sensor~chemical_condition, labeller = label_parsed,scales="free_y")
															            #+
															             


															save_plot(filename=paste0(graphTitle,"dualRecording_new_dFF_v3"),plot=tracePlot, plot_width=plot_dim*2.25,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 300)
															

															   
					#                             if(i == 1){
					# #### save to global env
					#                              traceB  <<- tracePlot

					#                             } else if(i == 2){
					#                              trace_APV_F <<- tracePlot
					                             	
                    # 							} else if(i == 3){

                    # 							}
				



					}
tracePlot 
}