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


path = "Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/"


# set working directory
dataDir = "Y:\\Sam/paper1_datasets/spont_v1/.csv/annotated_csv_activity_vs_marker"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)

singleAP_mainDir <- "Y:\\Sam/paper1_datasets/spont_v1/Figures/spont_figs_activity_vs_marker_v4"
mainDir = singleAP_mainDir







#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


full_file_list = list.files(pattern='.csv')

BasalRP = full_file_list[which(str_detect(full_file_list, pattern = "spont")==TRUE)]
#PairedPulse = full_file_list[which(str_detect(full_file_list, pattern = "singleAP|PP")==TRUE)]
#slowHz = full_file_list[which(str_detect(full_file_list, pattern = "1Hz|2Hz|5Hz")==TRUE)]
#fastHz = full_file_list[which(str_detect(full_file_list, pattern = "10Hz|20Hz")==TRUE)]
#spont = full_file_list[which(str_detect(full_file_list, pattern = "spont")==TRUE)]

groupings_list = list(BasalRP)
existing_groupings = groupings_list[lengths(groupings_list) > 0]



## run analysis routine across different dataset styles
for (item in existing_groupings) {

				if(exists("dataDir") ) {
					source(paste0(path,"peakAnalysis/set_dataDirectory.R") )
					} else {
						print("dataDir got cleared somewhere, expect an error.")
					}
				


	filenames = full_file_list[full_file_list %in% item]

	print(paste0("The following files will be analyzed : ", filenames) )

	
}








					print("Running stereotyped module for peak annotation, visualization of traces, and averaged trace outputs. This module is currently functional for datasets including BasalRP, PairedPulse, and 25 AP train stimulus paradigms.")

					source(paste0(path,"peakAnalysis/annotated_csv_reader.R") )
					

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=data
					 rm(data)
					 some_conditions = c("sensor","marker","segmentation","dish", "plate", "region","Ca", "protocol","exposeNum","ROINumber")
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "ROI_key" 
					 source(paste0(path,"peakAnalysis/get_mutate_condition.R") )
					 rm(df,some_conditions,other_conditions,.varname1)

 					
##### BEGIN ESTABLISHING DATASET CHARACTERISTICS, AUGMENTING DATAFRAME WITH PEAK IDENTIFIERS e.g. sync vs. async #####

					#should add the functionality  - if "stimParadigm" exists, then run getStimList
					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 		 
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					ROI_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum", "ROINumber")
					window_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum","ROINumber","windowedPeakID") #will be sunset in next pipeline
					num_APs = 25 # number of action potentials in a train protocol 
					s88x_offset = 0.05 #50 ms offset 



					 df_mutated = df %>% group_by(ROI_key) %>% mutate(record_time = max(absoluteTime) ) %>% dplyr::filter(windowedPeakID != "NotPeak")
					rm(df)
toc()


					#this could be set outside the script
					df = df_mutated
          some_conditions = c("sensor","dish","plate", "region","ROINumber" )
					  other_conditions=NULL
					 .varname1 = "trackROI_key" 
					 source(paste0(path,"peakAnalysis/get_mutate_condition.R") )
					 rm(df,some_conditions,other_conditions,.varname1, df_stimEpoch)




					df = df_mutated
					#rm(df_mutated)
					 some_conditions = c("dish","plate", "region","exposeNum")
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "vid_key" 
					 source(paste0(path,"peakAnalysis/get_mutate_condition.R") )
					 rm(df,some_conditions,other_conditions,.varname1)

					df = df_mutated 
					rm(df_mutated)
					#df = df %>% group_by(ROI_key,windowedPeakID) %>% mutate(normTime = absoluteTime - first(absoluteTime))
  					groupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 	
  					source(paste0(path,"peakAnalysis/expFitter_spont.R") )
  
  					df = fitPeaks
  					rm(fitPeaks)
					

					

					


  					# keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
  					# num_APs = 25 # number of action potentials in a train protocol 
  					# s88x_offset = 0.05 #50 ms offset 
  					# df = df_mutated
  					
  					# source(paste0(path,"peakAnalysis/getStimList_v2.R"))
  					# rm(keys, num_APs, s88x_offset)
  					



####### in this section of code, we will generate the necessary plots to compare 


					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber",'peakID')
					subgroupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","vid_key") 		 
					drops = c('data', '.resid','fit',"(weights)",".fitted")
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					
						#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_peakStats_spont_phys_v1.R")
					source(paste0(path,"peakAnalysis/get_peakStats_spont_phys_v1.R") )
          rm(drops)
  

					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					
					
					source(paste0(path,"peakAnalysis/get_plot_peakStats.R") )
					file_prefix = "activity_vs_marker_stats_v2_"
			#		source(paste0(path,"peakAnalysis/get_written_output_csv.R") )


					#empty_plot_area = NULL
  				#map_color = NULL
  				#map_C = NULL 
	 				margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
					gs = list(map_color, #2
								 map_B, #3
								 ROI_per_vid, 		#4
					       avg_traces,	#6
					       amplitude_histo,	#7
					       tau_histo,	#8
					       thalf_histo)	#9
					

					 hlay <- rbind(c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA),
					               c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA),
					               c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA),
					               c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA),
					               c(1,1,1,1,1,2,2,2,2,2,NA,NA,NA),
					               
					               c(3,3,3,3,4,4,4,4,4,4,4,4,4),
					               c(3,3,3,3,4,4,4,4,4,4,4,4,4),
					               c(3,3,3,3,4,4,4,4,4,4,4,4,4),
					               c(3,3,3,3,4,4,4,4,4,4,4,4,4),
					               
					               c(5,5,5,5,6,6,6,6,7,7,7,7,NA),
					               c(5,5,5,5,6,6,6,6,7,7,7,7,NA),
					               c(5,5,5,5,6,6,6,6,7,7,7,7,NA),
					               c(5,5,5,5,6,6,6,6,7,7,7,7,NA),
					               c(5,5,5,5,6,6,6,6,7,7,7,7,NA)
					               
					               )
				
				

    			scalefactor=0.75
					 	
					SIFigure1 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
				        #ggsave(filename="FigureS2_activityVSmarker_O-phys_paper_v1.emf",plot=SIFigure1,dpi=900, units="in",width=scalefactor*24,height=scalefactor*26, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="FigureS2_activityVSmarker_O-phys_paper_v1.jpeg",plot=SIFigure1,dpi=600, units="in",width=scalefactor*24,height=scalefactor*26, device = "jpeg")
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

