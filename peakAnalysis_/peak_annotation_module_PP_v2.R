#peak_annotation_module.R

#Written by Samuel T Barlow
# current as of 10/4/22

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

path = "Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/"





# set working directory
dataDir = "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/.csv/annotated_csv"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)

PP_mainDir <- "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/Figures/PP_v7_final"
mainDir = PP_mainDir



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
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
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




					#this could be set outside the script
					
					
					# ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
					# n=100 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
					# hide_titles=FALSE

					# source(paste0(path, "peakAnalysis/tracePlotter.R")) 
					# rm(ROI_IDs)
					



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
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "trackROI_key" 
					 source(paste0(path, "peakAnalysis/get_mutate_condition.R"))
					 rm(df,some_conditions,other_conditions,.varname1, df_stimEpoch)




					





					df = df_mutated
					rm(df_mutated)
					 some_conditions = c("dish","plate", "region" )
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
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
		 			levels = c(#'0pt25Ca',
		 						'0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c(#'#FDE725FF',
										#'#7AD151FF',
										'#22A884FF','#414487FF','#440154FF')
		



					ROIs_to_save = c("GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0005",
					                 "GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009")#,
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0003",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0004",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0005",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0006",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0007",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0008",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0010",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0011",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0012",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0013",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0014",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0015",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0016",
					                 #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0017")
					                 
					  
					  
					  
					  
					                #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0001",
					               #  "GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0002",
					                # "GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0003",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0004",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0005",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0006",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0007",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0008",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0009",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0010",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0011",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0012",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0013",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0014",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0015",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0016",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0017",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0018",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0019",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0020",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0021",
					                 ##"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0022",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0023",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0024",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0025",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0026",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0027",
					                 #"GluSnFR3-SyPhy-marker-dish05-plate10-region02-ROI0028")
					                 
					  
					  
					  #"GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0030",
										#"GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0033")#,
										#"GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0047")#, 
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0002",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0003",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0004",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0005",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0006",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0007",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0008",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0009",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0010",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0011",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0012",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0013",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0014",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0015",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0016",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0017",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0018",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0019",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0020",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0021",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0022",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0023",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0024",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0025",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0026",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0027",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0028",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0029",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0030",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0031",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0032",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0033",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0034",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0035",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0036",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0037",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0038",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0039",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0040",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0041",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0042",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0043",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0044",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0045",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0046",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0047",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0048",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0049",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0050",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0051",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0052",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0053",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0054",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0055",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0006",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0007",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0008",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0010",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0013",
										# #"GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0014",
										# "GluSnFR3-SyPhy-marker-dish06-plate01-region02-ROI0056")
					keep_PP = c("PP60","PP75","PP100", "PP150", "PP500")
					#tag_var_list = c('A','B','C','D')
					



  
  
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
  					#### I've changed the tag_var_list for the BRAIN poster!
            			# empty_plot_area = NULL
  					# empty_too = NULL
  					# ROI_plots = NULL
  					# subtraction = NULL
  					scalefactor=0.75
  tic()
  					print("try PPR plot")		
  					source(paste0(path, "peakAnalysis/get_PPR_plotter_v4_Figure4.R"))
  					
					#source(paste0(path, "peakAnalysis/get_PPR_plotter_figure5.R")
					#print("tried PPR plot")
		
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
				        ggsave(filename="Figure5_v7_dish6_plate01_region01_PPRphysiology_from_tracking_O-phys_paper.emf",plot=Figure5,dpi=900, units="in",width=30*scalefactor,height=42*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="Figure5_v7_dish6_plate01_region01__PPRphysiology_from_tracking_O-phys_paper.jpeg",plot=Figure5,dpi=600, units="in",width=30*scalefactor,height=42*scalefactor, device = "jpeg")


	 margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
				avg_PPR = NULL

					 gs = list(  map_output_1_PP100_PPR_bias_numeric)#,
					              #tmp_tracePlot1, #7
					             #tmp_PPRplot1, #8
					             #tmp_tracePlot2,#9
					             #tmp_PPRplot2)#, #8
					             

					 hlay <- rbind(
					               c(1,1,1,1,1,1),
					               c(1,1,1,1,1,1),
					               c(1,1,1,1,1,1),
					               c(1,1,1,1,1,1),
					               c(1,1,1,1,1,1),
					               c(1,1,1,1,1,1) 
					               											#12x25
					               )
				
				

    				
					Panel2_Figure4_partA_BRAIN = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
				        ggsave(filename="Panel2_Figure4_partA_BRAIN.emf",plot=Panel2_Figure4_partA_BRAIN,dpi=900, units="in",width=12*scalefactor,height=12*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="Panel2_Figure4_partA_BRAIN.jpeg",plot=Panel2_Figure4_partA_BRAIN,dpi=1200, units="in",width=12*scalefactor,height=12*scalefactor, device = "jpeg")


					 gs = list(  #map_output_1_PP100_PPR_bias_numeric)#,
					              tmp_tracePlot1, #7
					             tmp_PPRplot1, #8
					             tmp_tracePlot2,#9
					             tmp_PPRplot2)#, #8
					             

					 hlay <- rbind(
					               c(1,1,1,1,1,1,1,1,2,2,2,2),
					               c(1,1,1,1,1,1,1,1,2,2,2,2),
					               c(1,1,1,1,1,1,1,1,2,2,2,2),
					               c(1,1,1,1,1,1,1,1,2,2,2,2),
					               
					               c(3,3,3,3,3,3,3,3,4,4,4,4),
					               c(3,3,3,3,3,3,3,3,4,4,4,4),
					               c(3,3,3,3,3,3,3,3,4,4,4,4), 
					               c(3,3,3,3,3,3,3,3,4,4,4,4) 
					               											#12x25
					               )
				
				

    				
					Panel2_Figure4_partB_BRAIN = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
				        ggsave(filename="Panel2_Figure4_partB_BRAIN.emf",plot=Panel2_Figure4_partB_BRAIN,dpi=900, units="in",width=24*scalefactor,height=16*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
				        ggsave(filename="Panel2_Figure4_partB_BRAIN.jpeg",plot=Panel2_Figure4_partB_BRAIN,dpi=600, units="in",width=24*scalefactor,height=16*scalefactor, device = "jpeg")










toc()





				    # margin = theme(plot.margin = unit(c(2,2,2,2), "cm"))
				    
				    #    select_grobs <- function(lay) {
				    #   id <- unique(c(t(lay))) 
				    #   id[!is.na(id)]
				    # } 
					
					

toc()

					print("Finished stereotyped annotation module.")














