##pseudo-code_for_activity_marker_comparison.R


# written by Samuel T Barlow
# 12.5.23











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





library(tictoc)
tic()

# set working directory
dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/annotated_csv_activity_vs_marker_v2"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)

singleAP_mainDir <- "Y:\\Sam/paper1_datasets/singleAP_v2/Figures/activity_vs_slice_marker_v5_stats"
mainDir = singleAP_mainDir







#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


full_file_list = list.files(pattern='.csv')

BasalRP = full_file_list[which(str_detect(full_file_list, pattern = "singleAP")==TRUE)]
#PairedPulse = full_file_list[which(str_detect(full_file_list, pattern = "singleAP|PP")==TRUE)]
#slowHz = full_file_list[which(str_detect(full_file_list, pattern = "1Hz|2Hz|5Hz")==TRUE)]
#fastHz = full_file_list[which(str_detect(full_file_list, pattern = "10Hz|20Hz")==TRUE)]
#spont = full_file_list[which(str_detect(full_file_list, pattern = "spont")==TRUE)]

groupings_list = list(BasalRP)
existing_groupings = groupings_list[lengths(groupings_list) > 0]



## run analysis routine across different dataset styles
for (item in existing_groupings) {

				if(exists("dataDir") ) {
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/set_dataDirectory.R")
					} else {
						print("dataDir got cleared somewhere, expect an error.")
					}
				


	filenames = full_file_list[full_file_list %in% item]

	print(paste0("The following files will be analyzed : ", filenames) )

	
}








					print("Running stereotyped module for peak annotation, visualization of traces, and averaged trace outputs. This module is currently functional for datasets including BasalRP, PairedPulse, and 25 AP train stimulus paradigms.")

					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/annotated_csv_reader.R")
					

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=data
					 rm(data)
					 some_conditions = c("sensor","marker","segmentation","dish", "plate", "region","Ca", "protocol","exposeNum","ROINumber")
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "ROI_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)

 					
##### BEGIN ESTABLISHING DATASET CHARACTERISTICS, AUGMENTING DATAFRAME WITH PEAK IDENTIFIERS e.g. sync vs. async #####

					#should add the functionality  - if "stimParadigm" exists, then run getStimList
					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 		 
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					ROI_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum", "ROINumber")
					window_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum","ROINumber","windowedPeakID") #will be sunset in next pipeline
					num_APs = 25 # number of action potentials in a train protocol 
					s88x_offset = 0.051 #50 ms offset 


					df = df_mutated




					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/getStimList_v2.R")
					
					#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_sync_async.R")
					rm(window_keys,ROI_keys,groupers, num_APs, s88x_offset)




					#this could be set outside the script
					
					
					# ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
					# n=100 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
					# hide_titles=FALSE

					# source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/tracePlotter.R") 
					# rm(ROI_IDs)
					



					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=df_mutated
					 rm(df_mutated)
					 some_conditions = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "stimKey" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)





					 df = df_mutated
					 rm(df_mutated)



					keysEpoch = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber",'stimEpoch')
					groupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber") 
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_stimEpochs.R") 

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=df_stimEpoch
					 rm(df_stimEpoch)
					 some_conditions = c("sensor","marker","segmentation","dish","plate", "region","ROINumber" )
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "trackROI_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1, df_stimEpoch)




					





					df = df_mutated
					rm(df_mutated)
					 some_conditions = c("segmentation","dish","plate", "region","exposeNum" )
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "vid_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)

					df = df_mutated 
					rm(df_mutated)
					groupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 	
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/expFitter.R")

					df = fitPeaks
					rm(fitPeaks)

					# drops = c('data', '.resid','fit',"(weights)",".fitted")

					# source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_peak_halfwidths_rise_decaytimes_v1.R")


				# 	ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
				# 	n=50 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
				# 	hide_titles=TRUE
			       # add_max = TRUE
			       # add_time_dots = TRUE
				# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/tracePlotter.R") 
				# 	rm(ROI_IDs,add_max,hide_titles,n)



####### in this section of code, we will generate the necessary plots to compare 


					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber",'peakID')
					subgroupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","vid_key") 		 
					drops = c('data', '.resid','fit',"(weights)",".fitted")
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_peakStats_activity_vs_marker.R")
          rm(drops)
  

					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					
					
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_plot_peakStats.R")
					file_prefix = "activity_vs_marker_stats_v2_"
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_written_output_csv.R")


					empty_plot_area = NULL#ggplot() + 
  					# 					geom_rect(aes(xmin = 1, xmax = 3, ymin = 10, ymax = 15), fill = "white", alpha = 0.4, color = "white") + 
  					# 					labs(tag="A")+
  					# 					theme_void()+
  					#map_color = NULL	
	 				margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
					gs = list(empty_plot_area, 	#1
								 map_color, #2
								 map_C, #3
								 #map_C, #4
								 #map_D, #5
								 #map_E, #6
					             ROI_per_vid, 		#4
					             peak_probability,		#5
					             avg_traces,	#6
					             amplitude_histo,	#7
					             tau_histo,	#8
					             dt_histo)	#9
					

					 hlay <- rbind(#c(1,1,NA,NA,NA,NA,3,3,3,3,4,4,4,4),
					               #c(1,1,NA,NA,NA,NA,3,3,3,3,4,4,4,4),
					               #c(1,1,2,2,2,2,3,3,3,3,4,4,4,4),
					               #c(1,1,2,2,2,2,3,3,3,3,4,4,4,4),
					               c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),
					               c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),
					               c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),
					               c(1,1,1,1,1,1,2,2,2,2,3,3,3,3),
					                    
					               c(4,4,4,4,5,5,5,5,6,6,6,6,6,6),
					               c(4,4,4,4,5,5,5,5,6,6,6,6,6,6),
					               c(4,4,4,4,5,5,5,5,6,6,6,6,6,6),
					               c(4,4,4,4,5,5,5,5,6,6,6,6,6,6),
					               c(4,4,4,4,5,5,5,5,6,6,6,6,6,6), #5x14
					               
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA),
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA),
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA),
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA),
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA), #6x14
					               c(7,7,7,7,8,8,8,8,9,9,9,9,NA,NA)
					               
					               
					               )
				
				

    			scalefactor=0.75
					 	
					SIFigure1 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
				        ggsave(filename="Figure2_activityVSmarker_O-phys_paper_v6.emf",plot=SIFigure1,dpi=900, units="in",width=scalefactor*28,height=scalefactor*30, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="Figure2_activityVSmarker_O-phys_paper_v6.jpeg",plot=SIFigure1,dpi=600, units="in",width=scalefactor*28,height=scalefactor*30, device = "jpeg")
				        #ggsave(filename="Figure2_activityVSmarker_O-phys_paper_v3.png",plot=SIFigure1,dpi=600, units="in",width=scalefactor*28,height=scalefactor*38, device = "png")
				        


					



#### first, read in the necessary datasets

##annotated_csv_reader.R

### Then, apply peak annotation module stuff


### mutate_condition.R
### stimList.R
### def_stim_vlines.R


### then, apply expFitter and plot outputs to be sure. 

### expfitter.R

### tracePlotter.R


### Finally, generate comparisons of activity vs marker-based segmentation



### data cleaning. remove ROIs with zero trace (e.g. drift correction ROIs with no data from marker-based segmentation)
### ROIs found per video. should increase with activity based segmentation with Ca2+ 
### ROIs with a positive signal 
### comparing amplitude, exp decay , and dt with different segmentation methods. Add in statistical measrues, should be no difference hopefully 

