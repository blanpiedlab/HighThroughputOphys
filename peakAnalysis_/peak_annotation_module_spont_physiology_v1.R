##peak_annotation_module_singleAP_physiology_v2.R

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
dataDir = "Y:\\Sam/paper1_datasets/spont_v1/.csv/annotated_csv"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)

singleAP_mainDir <- "Y:\\Sam/paper1_datasets/spont_v1/Figures/spont_figs_v3_SI"
mainDir = singleAP_mainDir







#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


full_file_list = list.files(pattern='.csv')

#BasalRP = full_file_list[which(str_detect(full_file_list, pattern = "singleAP")==TRUE)]
#PairedPulse = full_file_list[which(str_detect(full_file_list, pattern = "singleAP|PP")==TRUE)]
#slowHz = full_file_list[which(str_detect(full_file_list, pattern = "1Hz|2Hz|5Hz")==TRUE)]
#fastHz = full_file_list[which(str_detect(full_file_list, pattern = "10Hz|20Hz")==TRUE)]
spont = full_file_list[which(str_detect(full_file_list, pattern = "spont")==TRUE)]

groupings_list = list(spont)
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





tic()


					print("Running stereotyped module for peak annotation, visualization of traces, and averaged trace outputs. This module is currently functional for datasets including BasalRP, PairedPulse, and 25 AP train stimulus paradigms.")

					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/annotated_csv_reader.R")
toc()			

					# #set up a script that adds mutated 2-condition secondAxis to the dataframe 
tic()
           df=data
					 rm(data)
toc()
tic()
					 some_conditions = c("sensor","marker","segmentation","dish", "plate", "region","Ca", "protocol","exposeNum","ROINumber")
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "ROI_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)

 					
##### BEGIN ESTABLISHING DATASET CHARACTERISTICS, AUGMENTING DATAFRAME WITH PEAK IDENTIFIERS e.g. sync vs. async #####

					
					df = df_mutated %>% group_by(ROI_key) %>% mutate(record_time = max(absoluteTime) ) %>% dplyr::filter(windowedPeakID != "NotPeak")
					rm(df_mutated)
toc()


					#this could be set outside the script
					
					
					 

					#df = df_mutated
					#rm(df_mutated)
					 some_conditions = c("dish","plate", "region","exposeNum")
					
					 #c("dish", "plate", "region", "replicate", "ROINumber")	
												
						
					 other_conditions=NULL
					 .varname1 = "vid_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)

					df = df_mutated 
					rm(df_mutated)
					#df = df %>% group_by(ROI_key,windowedPeakID) %>% mutate(normTime = absoluteTime - first(absoluteTime))
  					groupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 	
  					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/expFitter_spont.R")
  
  					df = fitPeaks
  					rm(fitPeaks)
  					
  				
												
  				 some_conditions = c("sensor","dish","plate", "region","Ca", "ROINumber" )
					 other_conditions=NULL
					 .varname1 = "trackROI_key" 
					 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_mutate_condition.R")
					 rm(df,some_conditions,other_conditions,.varname1)

  
  					
  					# drops = c('data', '.resid','fit',"(weights)",".fitted")
  					# source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_light_df.R") 
  					# orig_df = df
  					# rm(df)
  					
  					
           			# df = light_df 
            		#rm(light_df)
                        
           	# groupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber","ROI_key",'peakID','windowedPeakID')
  			# 		drops = c('data', '.resid','fit',"(weights)",".fitted")
  					
            # 		source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_peak_halfwidths_rise_decaytimes_v2_spont.R")
  

  			# 		ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
  			# 		 n=25 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
  			# 		 hide_titles=FALSE
            #  			 add_max = TRUE
  			# 		 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/tracePlotter_spont_v2.R") 
  			# 		 rm(ROI_IDs)

  			# 		 ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
  			# 		 n=25 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
  			# 		 hide_titles=FALSE
            #  add_max = TRUE
            #  add_time_dots = TRUE
  			# 		 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/peakPlotter_v2_spont.R") 
  			# 		 rm(ROI_IDs)




###### in this section of code, we will generate the necessary plots to compare 
			    	
			    	df= df_mutated
			    	groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key","trackROI_key","ROINumber",'peakID')
					subgroupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","vid_key","trackROI_key") 		 
					drops = c('data', '.resid','fit',"(weights)",".fitted")
					#keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_peakStats_spont_phys_v1.R")
          #rm(drops)
  				#rm(df_mutated)








  				
  

					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					
					
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_plot_peakStats_spont_v1.R")
					file_prefix = "spont_phys_v1_"
					source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/get_written_output_csv.R")


# # 					empty_plot_area = NULL
# 	 				margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
# 					gs = list(empty_plot_area, 	#1
# 					             ROI_per_vid, 		#2
# 					             peak_probability,		#3
# 					             avg_traces,	#4
# 					             amplitude_histo,	#5
# 					             tau_histo,	#6
# 					             dt_histo)	#7
					

# 					 hlay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1),
# 					               c(1,1,1,1,1,1,1,1,1,1,1,1),
# 					               c(1,1,1,1,1,1,1,1,1,1,1,1),
					               
# 					               c(2,2,2,3,3,3,4,4,4,4,4,4),
# 					               c(2,2,2,3,3,3,4,4,4,4,4,4),
# 					               c(2,2,2,3,3,3,4,4,4,4,4,4),
# 					               c(2,2,2,3,3,3,4,4,4,4,4,4),
					               
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7),
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7),
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7),
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7),
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7),
# 					               c(5,5,5,5,6,6,6,6,7,7,7,7)
					               
					               
# 					               )
				
				

    				
# 					SIFigure1 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
# 				        ggsave(filename="SIFigure1_activityVSmarker_O-phys_paper.emf",plot=SIFigure1,dpi=900, units="in",width=24,height=30, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
# 				        ggsave(filename="SIFigure1_activityVSmarker_O-phys_paper.jpeg",plot=SIFigure1,dpi=600, units="in",width=24,height=30, device = "jpeg")



					

