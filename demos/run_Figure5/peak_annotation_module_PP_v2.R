#peak_annotation_module.R

#Written by Samuel T Barlow
# current as of 2/19/2025

# At present, this module goes through most of the typical paces of an analysis session.
# From annotated, per stimulus paradigm, "master" .csv files with peaks identified, we can:
# 1) perform an exponential decay fit of every identified peak and ensure that each fit passes a quality check (monotonically decreasing)
# 2) apply additional identifiers to the data.frame, including "ROI_keys" (for sorting procedures that depend on unique ROIs, e.g. what was the cumulative dF/F at each ROI according to stimulus paradigm?),
#																"stimKeys"(for characterizing stimulus trains on a per video basis),
#																"epoch_keys" (for characterizing what is going on in each stimulus window)
# 3) Output a subset of all traces along with characterized peaks. Currently, this is set to 50. It's possible that this should be set by the user outside of the script shell.
# 4) mutate a specific sorting column(s) to allow for an interaction. Here, I wanted to sort data based on neuron segment (axon or dendrite) AND Ca2+ condition (0.5mMCa,1mMCa...) for visualizing how traces change with each defined condition. This is flexible. 
# 5) Averaged traces. averageTraces.R takes the data.frame, chops off the first few frames to get us right up next to the stimulus, and then averages all the existing traces together while plotting the unique traces in the background (by ROI_keys). 
#	 This is handy because you can see the macro-level behavior of the neurons under different bath/imaging/protocol conditons. In the background, I've added a color-coded ribbon that tells you which stimulus epoch we're in - this is handy because we can see 
#	 the behavior of the stimEpochs assignment. It works! 


# At this point, the module has been explicitly tested across the following stimulus paradigms: BasalRP, PP100, PP500, PP1000, 1Hz, 2Hz, 5Hz. I anticipate that it will work without any errors on PP50 and PP2500. 
# This script will not work on 10Hz or 20Hz (more massaging needed there.).



# The current final output of this script is the data.frame, fitsLabeled_timeClass_stimKeyed_mutated. 

# Currently omitted from this script is the half-width determination for every peak, and individual peak plotting, which go hand in hand. With the present dataset, half-width determinations are problematic because of the TTL_start effect I am seeing:
# briefly, when Fusion fires the TTL using the Spare Trigger, the frame at which this occurs takes 10 times longer. This is introducing imaging artifacts into the first peak in a stimulus paradigm, especially apparent in BasalRP paradigms. 
# So I will not be re-implementing half_widths.R until I have fixed my acquisitions.   







#tic()

##### if this is your first session in R, please uncomment the install.packages() list

#list_of_packages = c('plyr','dplyr','ggplot2','zoo','ggthemes','tidyr','reshape2','stringr','readr','purrr','extrafont','tictoc','broom','gridExtra','grid','imager','readr','EBImage','bioimagetools')

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


## Specify the path to the script directory here. 

path = "Z:\\R Scripts/Github_repo/"





# set working directory
dataDir = "Z:\\R Scripts/Github_repo/demos/run_Figure5/annotated_csv"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)


###set your figure directory

PP_mainDir <- "Z:\\R Scripts/Github_repo/demos/run_Figure5/Figures_dir"
mainDir = PP_mainDir
png_dir = "Z:\\R Scripts/Github_repo/demos/run_Figure5/zstackOutputs"




tic()


#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


full_file_list = list.files(pattern='.csv')

BasalRP = full_file_list[which(str_detect(full_file_list, pattern = "BasalRP")==TRUE)]
PairedPulse = full_file_list[which(str_detect(full_file_list, pattern = "singleAP|PP")==TRUE)]
slowHz = full_file_list[which(str_detect(full_file_list, pattern = "1Hz|2Hz|5Hz")==TRUE)]
fastHz = full_file_list[which(str_detect(full_file_list, pattern = "10Hz|20Hz")==TRUE)]
spont = full_file_list[which(str_detect(full_file_list, pattern = "spont")==TRUE)]

groupings_list = list(BasalRP,PairedPulse,slowHz,fastHz,spont)
existing_groupings = groupings_list[lengths(groupings_list) > 0]



## run analysis routine across different dataset styles
for (item in existing_groupings) {

				if(exists("dataDir") ) {
					source(paste0(path, "peakAnalysis/set_dataDirectory.R"))
					} else {
						print("dataDir got cleared somewhere, expect an error.")
					}
				


	filenames = full_file_list[full_file_list %in% item]

	print(paste0("The following files will be analyzed : ", filenames) )

	
}

#print("The analysis loop should be done.")	








					print("Running stereotyped module for peak annotation, visualization of traces, and averaged trace outputs. This module is currently functional for datasets including BasalRP, PairedPulse, and 25 AP train stimulus paradigms.")

					source(paste0(path, "peakAnalysis/annotated_csv_reader.R"))
					

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=data
					 rm(data)
					 some_conditions = c("sensor","marker","segmentation","TTL_start","dish", "plate", "region","Ca", "protocol","exposeNum","ROINumber")
					 							
						
					 other_conditions=NULL
					 .varname1 = "ROI_key" 
					 source(paste0(path, "peakAnalysis/get_mutate_condition.R"))
					 rm(df,some_conditions,other_conditions,.varname1)

 					
##### BEGIN ESTABLISHING DATASET CHARACTERISTICS, AUGMENTING DATAFRAME WITH PEAK IDENTIFIERS e.g. sync vs. async #####

					#should add the functionality  - if "stimParadigm" exists, then run getStimList
					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 		 
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					ROI_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum", "ROINumber")
					window_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum","ROINumber","windowedPeakID") #will be sunset in next pipeline
					num_APs = 25 # number of action potentials in a train protocol 
					s88x_offset = 0.050 #50 ms offset 


					df = df_mutated




					source(paste0(path, "peakAnalysis/getStimList_v2.R"))
					
					rm(window_keys,ROI_keys,groupers)




					



					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=df_mutated
					 rm(df_mutated)
					 some_conditions = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "stimKey" 
					 source(paste0(path, "peakAnalysis/get_mutate_condition.R"))
					 rm(df,some_conditions,other_conditions,.varname1)





					 df = df_mutated
					 rm(df_mutated)



					keysEpoch = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber",'stimEpoch')
					groupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber") 
					source(paste0(path, "peakAnalysis/get_stimEpochs.R") )

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
					 df=df_stimEpoch
					 rm(df_mutated)
					 some_conditions = c("sensor","marker","segmentation","dish","plate", "region","ROINumber" )
					
												
						
					 other_conditions=NULL
					 .varname1 = "trackROI_key" 
					 source(paste0(path, "peakAnalysis/get_mutate_condition.R"))
					 rm(df,some_conditions,other_conditions,.varname1, df_stimEpoch)




					





					df = df_mutated
					rm(df_mutated)
					 some_conditions = c("dish","plate", "region" )
					
						
					 other_conditions=NULL
					 .varname1 = "vid_key" 
					 source(paste0(path, "peakAnalysis/get_mutate_condition.R"))
					 rm(df,some_conditions,other_conditions,.varname1)





					df = df_mutated
					
					rm(df_mutated)
					groupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca","Ca_mM", "protocol",  "exposeNum","ROINumber","ROI_key","trackROI_key","vid_key")
					stimEpoch_groupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca","Ca_mM", "protocol",  "exposeNum","ROINumber","ROI_key","stimEpoch","trackROI_key","vid_key")
					keys = c('trackROI_key')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#22A884FF','#414487FF','#440154FF')
		


          ## this statement can be changed to select different ROIs for display. Try editing the ROI ID number!! ROI0005 and ROI0009 are the ones used in the manuscript. 
					ROIs_to_save = c("GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0005",
					                 "GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009")
					                 
					  
					  
					keep_PP = c("PP60","PP75","PP100", "PP150", "PP500")
					



  
  
           # 			Figure 5 implementation
  					empty_plot_area = NULL
  					avg_peaks = NULL
  					avg_PPR = NULL
  					violin_PPR = NULL
  					spat_color = NULL
  					spat_overlay = NULL
  					spat_facet = NULL
  					ROI_A_facet = NULL
  					ROI_B_facet = NULL
  
  						n=2             #1....2...3...4....5..6...7...8...9..10..11..12..13..14..15..16..17..18
  					tag_var_list = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R")
  					scalefactor=0.75
  tic()
  					print("try PPR plot")		
  					source(paste0(path, "peakAnalysis/get_PPR_plotter_v4.R"))
  					
			
				 margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
				avg_PPR = NULL

					 gs = list(  avg_peaks, 		#2
					             subtract_example,
					             violin_PPR,
					             PPR_bias,		#3
					             map_output_1_PP100_PPR_bias_numeric,
					             map_output_1_PP60_mean_PPR_perROI,	#4
					           	#5
					           map_output_1_PP100_mean_PPR_perROI,	#6
					             tmp_tracePlot1, #7
					             tmp_PPRplot1, #8
					             tmp_tracePlot2,#9
					             tmp_PPRplot2)#, #8
					             

					 hlay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), #4x14
					 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
					 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
					 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
					 			   
					 			   c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4),
					               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4),
					               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4), #3x13
					               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4), #3x13
					               
					               c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7),
					               c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7), #4x14
					               c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7),
					               c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7),
					               c(5,5,5,5,5,6,6,6,6,6,7,7,7,7,7),

					               c(8,8,8,8,8,8,8,8,8,8,8,9,9,9,9),
					               c(8,8,8,8,8,8,8,8,8,8,8,9,9,9,9), #4x14
					               c(8,8,8,8,8,8,8,8,8,8,8,9,9,9,9),
					               c(8,8,8,8,8,8,8,8,8,8,8,9,9,9,9),
					               
					               c(10,10,10,10,10,10,10,10,10,10,10,11,11,11,11),
					               c(10,10,10,10,10,10,10,10,10,10,10,11,11,11,11),
					               c(10,10,10,10,10,10,10,10,10,10,10,11,11,11,11),
					               c(10,10,10,10,10,10,10,10,10,10,10,11,11,11,11) #4x14
					                
					               											#12x25
					               )
				
				

    				
					Figure5 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
				        ggsave(filename="Figure5_O-phys_paper.emf",plot=Figure5,dpi=900, units="in",width=30*scalefactor,height=42*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="Figure5_O-phys_paper.jpeg",plot=Figure5,dpi=600, units="in",width=30*scalefactor,height=42*scalefactor, device = "jpeg")



					
					

toc()

					print("Finished stereotyped annotation module.")














