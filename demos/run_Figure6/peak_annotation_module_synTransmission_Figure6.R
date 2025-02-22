#peak_annotation_module_mini.R

#Written by Samuel T Barlow
# current as of 2/19/25




##### if this is your first session in R, please uncomment the install.packages() list

#list_of_packages = c('plyr','dplyr','ggplot2','zoo','ggthemes','tidyr','reshape2','stringr','readr','purrr','extrafont','tictoc','broom','gridExtra','grid','imager','readr','EBImage','bioimagetools','rstatix','nls.multstart','nlstools','devtools')

#install.packages(list_of_packages)


#necessary libraries
library(plyr) ###File manipulation
library(dplyr) ###File manipulation
library(ggplot2) ###plotting
library(zoo)  ###not used in this, but mathematical functions
library(ggthemes) ###other ggplot2 graph themes
library(tidyr)
library(reshape2)
library(stringr)
library(readr)
library(ggpubr)
library(purrr)
library(extrafont)
library(tictoc)
library(broom)
library(gridExtra)
library("RColorBrewer")
library(rstatix)

scalefactor=1




#set path for R scripts
path = "YourDrive:\\HighThroughputOphys-main/"


# set directory for annotated_csv
# Please find the data at the bottom of the README
dataDir = "YourDrive:\\HighThroughputOphys-main/annotated_csv_Figure6_analysis"
setwd(dataDir)
filenames = list.files(pattern=".csv")


#set directory for figures
mainDir <- "YourDrive:\\HighThroughputOphys-main/demos/run_Figure6/Figures_dir"

#set path to find the .png files for spatial mapping
png_path = "YourDrive:\\HighThroughputOphys-main/demos/run_Figure6/png_dir"

tic()



					print("Running stereotyped module for peak annotation, visualization of traces collected from spontaneous dual-color recordings.")

					source(paste0(path,"peakAnalysis/annotated_csv_reader.R"))

					df=data
					rm(data)
					some_conditions = c('sensor',"dish", "plate", "region","chemical_condition", "ROINumber")
					#c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					other_conditions=NULL
					.varname1 = "conditionROI_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)

					df = df_mutated
					rm(df_mutated)
					some_conditions = c("sensor","dish","plate", "region","ROINumber")	
					other_conditions=NULL
					.varname1 = "sensor_trackROI_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)

					
					df = df_mutated
					rm(df_mutated)
					some_conditions = c("dish","plate", "region","ROINumber")	
					other_conditions=NULL
					.varname1 = "trackROI_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)


					df = df_mutated
					rm(df_mutated)
					some_conditions = c("dish","plate", "region" )	
					other_conditions=NULL
					.varname1 = "vid_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)

					

          
					
					plotBy = c("chemical_condition")
					levels = c("ctrl", "APV")
					color_override = c('grey21', 'blue')
					my_comparisons = NULL
					
					vid_keys_to_keep = c("dish10-plate04-region02",
					                     "dish7-plate10-region04",
					                     "dish7-plate08-region03",
					                     "dish7-plate08-region04",
					                     "dish7-plate07-region04",
					                     "dish7-plate07-region03",
					                     "dish8-plate08-region07",
					                     "dish8-plate08-region05",
					                     "dish8-plate08-region03",
					                     "dish8-plate07-region05")
					
					filtered_df_traces = df_mutated 
					variance_filter<- filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, 
					                                                       alpha_str == "notPeak",
					                                                      chemical_condition == "ctrl",
					                                                      sensor == "JF646") %>% 
					  group_by(trackROI_key) %>%
					  summarise(sd_trace = sd(dFF,na.rm=TRUE)) %>%
					  ungroup() #%>%
					  #dplyr::filter(sd_trace <= quantile(sd_trace, 0.85, na.rm=TRUE))
					filtered_df_traces<- filtered_df_traces %>% dplyr::filter(trackROI_key %in% unique(variance_filter$trackROI_key) )
					
					vid_keys_to_keep = c("dish10-plate04-region02",
					                     "dish7-plate10-region04",
					                     "dish7-plate08-region03",
					                     "dish7-plate08-region04",
					                     "dish7-plate07-region04",
					                     "dish7-plate07-region03",
					                     "dish8-plate08-region07",
					                     "dish8-plate08-region05",
					                     "dish8-plate08-region03",
					                     "dish8-plate07-region05")
					ROIs_to_save = c(#"dish7-plate07-region04-ROI0020", 
					                 "dish7-plate10-region04-ROI0010","dish8-plate08-region02-ROI0003")
					source(paste0(path,"peakAnalysis/ROIs_to_remove_synTrans_v1.R"))
        			


					# small_df_traces = filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked) #%>% 
        			#                                       #dplyr::filter(chemical_condition == "ctrl") #%>%
        			#                                       #dplyr::filter(trackROI_key == "dish7-plate07-region04-ROI0019")
					



        			# #### display all traces in dataset
					# tic()
					
					# peak_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition','peakID')
					# groupers = c('chemical_condition',
					#              'trackROI_key')
					# positive_ROI_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition')
					# plotBy = c("chemical_condition")
					# levels = c("ctrl", "APV")
					# color_override = c('grey21', 'blue')
					# my_comparisons = NULL
					# ROIs = unique(small_df_traces$trackROI_key)
					# n = 1
					# ROIs_to_sample = ceiling( length(ROIs) /n ) 
					# ROIs_to_save<- ROIs#sample(ROIs,ROIs_to_sample)
					# tmp_tag = NULL
					# disp_trans = TRUE
					# source(paste0(path,"peakAnalysis/get_single_ROI_display_v2.R"))
					# toc()
					#traceB<<-dual_record_display


        			small_df_traces = filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked) %>% 
        			                                     #dplyr::filter(vid_key == "dish8-plate08-region07")
        			                                      dplyr::filter(chemical_condition == "ctrl") %>%
        			                                      dplyr::filter(trackROI_key == "dish7-plate07-region04-ROI0019")
					

					
					
					peak_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition','peakID')
					groupers = c('chemical_condition',
					             'trackROI_key')
					positive_ROI_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition')
					plotBy = c("chemical_condition")
					levels = c("ctrl", "APV")
					color_override = c('grey21', 'blue')
					my_comparisons = NULL
					ROIs = unique(small_df_traces$trackROI_key)
					n = 1
					ROIs_to_sample = ceiling( length(ROIs) /n ) 
					ROIs_to_save<- ROIs#sample(ROIs,ROIs_to_sample)
					tmp_tag = c("B","C")
					disp_trans = FALSE
					source(paste0(path,"peakAnalysis/get_single_ROI_display_v2.R"))
					
					traceB<<-dual_record_display[[1]]
					scatter_amplitude_corr_replacement_C <<- dual_record_display[[2]]

					small_df_traces = filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked) %>% 
					  #dplyr::filter(chemical_condition == "ctrl") %>%
					  dplyr::filter(trackROI_key == "dish8-plate07-region05-ROI0020")
					
					



        			small_df_traces = filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked) %>% 
        			                                     #dplyr::filter(vid_key == "dish8-plate08-region07")
        			                                      #dplyr::filter(chemical_condition != "ctrl") %>%
        			                                      dplyr::filter(trackROI_key == "dish7-plate07-region04-ROI0019")
					

					
					
					peak_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition','peakID')
					groupers = c('chemical_condition',
					             'trackROI_key')
					positive_ROI_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition')
					plotBy = c("chemical_condition")
					levels = c("ctrl", "APV")
					color_override = c('grey21', 'blue')
					my_comparisons = NULL
					ROIs = unique(small_df_traces$trackROI_key)
					n = 1
					ROIs_to_sample = ceiling( length(ROIs) /n ) 
					ROIs_to_save<- ROIs#sample(ROIs,ROIs_to_sample)
					tmp_tag = c("","")
					disp_trans = TRUE
					source(paste0(path,"peakAnalysis/get_single_ROI_display_v2.R"))
					
					trace_SI_A<<-dual_record_display[[1]]
					scatter_amplitude_corr_replacement_SI_B <<- dual_record_display[[2]]



					margin = theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
				
			 
					gs = list(trace_SI_A#, 	#1
					          #scatter_amplitude_corr_replacement_SI_B #2
					          
					          )

					hlay<- rbind(
								
								c(1,1,1,1,1,1,1,1,1,1,1),#,2,2,2,2,2,2),
								c(1,1,1,1,1,1,1,1,1,1,1),#2,2,2,2,2,2),
								c(1,1,1,1,1,1,1,1,1,1,1),#2,2,2,2,2,2),
								c(1,1,1,1,1,1,1,1,1,1,1)#2,2,2,2,2,2)
								
								)

									#15x10
          scalefactor=0.75
    				library(devEMF)
					Figure6 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
					ggsave(filename="SI_Figure1_APV_activation.emf",plot=Figure6,dpi=900, units="in",width=22*scalefactor,height=8*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
					ggsave(filename="SI_Figure1_APV_activation.jpeg",plot=Figure6,dpi=600, units="in",width=22*scalefactor,height=8*scalefactor, device = "jpeg")
					




					
					
					peak_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition','peakID')
					groupers = c('chemical_condition',
					             'trackROI_key')
					positive_ROI_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition')
					plotBy = c("chemical_condition")
					levels = c("ctrl", "APV")
					color_override = c('grey21', 'blue')
					my_comparisons = NULL
					ROIs = unique(small_df_traces$trackROI_key)
					n = 1
					ROIs_to_sample = ceiling( length(ROIs) /n ) 
					ROIs_to_save<- ROIs#sample(ROIs,ROIs_to_sample)
					tmp_tag = c("C","D")
					disp_trans = TRUE
					source(paste0(path,"peakAnalysis/get_single_ROI_display_v2.R"))
					
					trace_vehicle_C<<-dual_record_display[[1]]
					scatter_amplitude_corr_replacement_D <<- dual_record_display[[2]]
					
					
					small_df_traces = filtered_df_traces %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked) %>% 
        			                                      #dplyr::filter(chemical_condition == "ctrl") %>%
        			                                      dplyr::filter(trackROI_key == "dish7-plate10-region04-ROI0016")
					

					
					
					peak_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition','peakID')
					groupers = c('chemical_condition',
					             'trackROI_key')
					positive_ROI_groupers = c('sensor','vid_key','trackROI_key','sensor_trackROI_key','conditionROI_key','chemical_condition')
					plotBy = c("chemical_condition")
					levels = c("ctrl", "APV")
					color_override = c('grey21', 'blue')
					my_comparisons = NULL
					ROIs = unique(small_df_traces$trackROI_key)
					n = 1
					ROIs_to_sample = ceiling( length(ROIs) /n ) 
					ROIs_to_save<- ROIs#sample(ROIs,ROIs_to_sample)
					tmp_tag = c("G","H")
					disp_trans = TRUE
					source(paste0(path,"peakAnalysis/get_single_ROI_display_v2.R"))
					
					trace_APV_G<<-dual_record_display[[1]]
					scatter_amplitude_corr_replacement_H <<- dual_record_display[[2]]
					


					
					


					vid_keys_to_keep = c("dish10-plate04-region02",
					                     "dish7-plate10-region04",
					                     "dish7-plate08-region03",
					                     "dish7-plate08-region04",
					                     "dish7-plate07-region04",
					                     "dish7-plate07-region03",
					                     "dish8-plate08-region07",
					                     "dish8-plate08-region05",
					                     "dish8-plate08-region03",
					                     "dish8-plate07-region05")
					ROIs_to_save = c(#"dish7-plate07-region04-ROI0020", 
					                 "dish7-plate10-region04-ROI0010","dish8-plate08-region02-ROI0003")
					source(paste0(path,"peakAnalysis/ROIs_to_remove_synTrans_v1.R"))
        
					df_traces = df_mutated %>% dplyr::filter(vid_key %in% vid_keys_to_keep, !trackROI_key %in% ROIs_to_remove_borked)


					
					ROIs_to_save = c("dish7-plate07-region04-ROI0004",
					                  "dish7-plate07-region04-ROI0020")
					
					source(paste0(path, "peakAnalysis/test_maxFinder_GluCa.R"))
					
					
					
					



#### FIGURE 6 IMPLEMENTATION



			
					
					margin = theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
				
	        empty_plot_area_A = NULL
	        empty_plot_area_J = NULL
	        scatter_Gluamplitude_dt_replacement_D = NULL #4
			histo_success_replacement_E = NULL #5
			peak_rate_replacement_H = NULL #8
			peak_transmission_replacement_I = NULL #9
					          
		 
					gs = list(empty_plot_area_A, 	#1
					          traceB, #2
					          scatter_amplitude_corr_replacement_C, #3
					          #trace_vehicle_D,
					          #scatter_amplitude_corr_replacement_E,
					          #trace_APV_F,
					          #scatter_amplitude_corr_replacement_G,
					          avg_traces_dualcolor,
					          amplitude_z_score_corr,
					          avg_traces_dualcolor_1,
					          avg_traces_dualcolor_2,
					          regionF,
					          tau_violin
					          #Glu_dFF_histo
					          #JF_dFF_histo

					          #empty_plot_area_J,
					          
					          )

					hlay<- rbind(
								c(1,1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3),
								
								c(2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3), #8x14

								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5),
								

								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9),
								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9),
								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9),
								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9),
								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9),
								c(8,8,8,8,8,8,6,6,6,7,7,7,9,9,9,9,9)		
								
								)

									#15x10
          scalefactor=0.75
    				library(devEMF)
					Figure6 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
					ggsave(filename="Figure6_v6-color_physiology_O-phys_paper.emf",plot=Figure6,dpi=900, units="in",width=34*scalefactor,height=30*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
					ggsave(filename="Figure6_v6-color_physiology_O-phys_paper.jpeg",plot=Figure6,dpi=600, units="in",width=34*scalefactor,height=30*scalefactor, device = "jpeg")
					



##### EXTENDED DATA FIGURE 6 MOCKUP


          empty_plot_area = NULL
								margin = theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
				
	     
					gs = list(empty_plot_area,
								region_vehicle_A, 	
					          region_vehicle_B, 
					          trace_vehicle_C,
					          scatter_amplitude_corr_replacement_D,					          
					          
					          empty_plot_area,

					          region_APV_E, 	
					          region_APV_F, 
					          trace_APV_G,
					          scatter_amplitude_corr_replacement_H
					          )

					hlay<- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),

								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								c(2,2,2,2,2,2,2,3,3,3,3,3,3,3),
								

								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5),
								c(4,4,4,4,4,4,4,4,5,5,5,5,5,5),
								
								c(6,6,6,6,6,6,6,6,6,6,6,6,6,6),

								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								c(7,7,7,7,7,7,7,8,8,8,8,8,8,8),
								
								c(9,9,9,9,9,9,9,9,10,10,10,10,10,10),
								c(9,9,9,9,9,9,9,9,10,10,10,10,10,10),
								c(9,9,9,9,9,9,9,9,10,10,10,10,10,10),
								c(9,9,9,9,9,9,9,9,10,10,10,10,10,10)
								
								
								
								)

									#15x10
          scalefactor=0.75
    				library(devEMF)
					Figure6 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
					ggsave(filename="ext_Figure6_v1-color_physiology_O-phys_paper.emf",plot=Figure6,dpi=900, units="in",width=28*scalefactor,height=48*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
					ggsave(filename="ext_Figure6_v1-color_physiology_O-phys_paper.jpeg",plot=Figure6,dpi=600, units="in",width=28*scalefactor,height=48*scalefactor, device = "jpeg")
					














toc()

					print("Finished stereotyped annotation module.")
					



