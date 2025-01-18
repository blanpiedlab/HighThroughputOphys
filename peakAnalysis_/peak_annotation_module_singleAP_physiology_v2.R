##peak_annotation_module_singleAP_physiology_v2.R

# written by Samuel T Barlow
# 12.5.23

#updated and commented 4.29.24

## peak_annotation_module_singleAP_physiology_v2.R is a collection of scripts for the generation of Figure 3 and Figure 4 in the manuscript: 
#"Dissecting  the function of single glutamatergic synapses with high-throughput optical physiology"











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






tic()

# set working directory
dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/annotated_csv_marker_only_v2"
#dataDir = "Y:\\Sam/paper1_datasets/singleAP_v2/.csv/subset_annotated_csv"

setwd(dataDir)

singleAP_mainDir <- "Y:\\Sam/paper1_datasets/singleAP_v2/Figures/singleAP_phys_v8_revisingFigures"
mainDir = singleAP_mainDir







#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


full_file_list = list.files(pattern='.csv')

BasalRP = full_file_list[which(str_detect(full_file_list, pattern = "singleAP")==TRUE)]

groupings_list = list(BasalRP)
existing_groupings = groupings_list[lengths(groupings_list) > 0]



## run analysis routine across different dataset styles
for (item in existing_groupings) {

				if(exists("dataDir") ) {
					source(paste0(path,"peakAnalysis/set_dataDirectory.R"))
					} else {
						print("dataDir got cleared somewhere, expect an error.")
					}
				


	filenames = full_file_list[full_file_list %in% item]

	print(paste0("The following files will be analyzed : ", filenames) )

	
}








					print("Running module for peak annotation, visualization of traces, and Figure generation. This module combines data from spontaneous and evoked recordings of iGluSnFR3 activity.")


## annotated_csv_reader.R is a script that reads in the annotated data from peakFinder.R

					source(paste0(path,"peakAnalysis/annotated_csv_reader.R"))
					
 					df=data
					rm(data)
					 


## mutate_condition.R is a basic function to append combined variables to the dataframe. mutate_condition.R accepts:
# .varname1, a column header for the new dataframe 
# some_conditions, a vector of strings that specify the column names to be pasted into a single variable.
# other_conditions, a deprecated variable that can be used to specify additional columns. 

# mutate_condition.R simply pastes together the data for each row from each column in some_conditions, separated by a "-".
# Example behavior. ROI_key, a new column with the following string: "GluSnFR3-SyPhy-mRuby3-marker-dish01-plate01-region01-1Ca-singleAP-repl01-ROI0001"	

					

					### adding ROI_key to dataframe
					some_conditions = c("sensor","marker","segmentation","dish", "plate", "region","Ca", "protocol","exposeNum","ROINumber")					 							
					other_conditions=NULL
					.varname1 = "ROI_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)





## getStimList_v2.R is a function that generates a list of key value pairs for the stimulus paradigm. In these experiments, the field stimulus was not directly measured. 
#  Rather, a TTL was delivered at a defined frame (Frame 195), and this knowledge was used to interpolate the stimulus position for every video.
#  NOTE: In our system, the TTL delivery coincided with single frame slow-down of the framerate (5 ms exposure time; normally 5.6-6.5 ms per frame, spiking to ~33 ms during TTL delivery). 
#        To avoid this imaging artifact, we manually added a time delay to the S88X Grass Stimulator, such that for every TTL it received from the Master-8, a 50 ms delay was added before delivering the stimulus to the neurons. 
#		 This value is encoded in s88x_offset. 

#  To interpolate the first stimulus, s88x_offset is added to the time stamp of the frame BEFORE the stimulus (i.e. Frame 194). 
#  If additional stimuli occur in the stimulus paradigm (e.g. paired pulse or train), the additional time information are extracted from the "protocol" string variable (e.g. "PP500" = interstimulus interval of 500 ms). 
#  For stimulus trains, an additional variable, num_APs is used to specify the length in stimuli of the train. These are distributed according to a frequency defined by the "protocol" string.     
 					

 					#### DEPRECATED; VERIFY CAN BE DELETED  ####

			 					#should add the functionality  - if "stimParadigm" exists, then run getStimList
								#groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 		
								#ROI_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum", "ROINumber")
								#window_keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum","ROINumber","windowedPeakID") #will be sunset in next pipeline
					 
					
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					num_APs = 25 # number of action potentials in a train protocol 
					s88x_offset = 0.05 #50 ms offset 
					df = df_mutated

					source(paste0(path,"peakAnalysis/getStimList_v2.R"))
					rm(keys, num_APs, s88x_offset)




					
					### adding stimKey to dataframe
					df=df_mutated
					rm(df_mutated)
					some_conditions = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					other_conditions=NULL
					.varname1 = "stimKey" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)




## get_stimEpochs.R is a function that defines stimulus epochs according to which protocol is being run. 
#  This is especially useful in paired pulse paradigms, where the trace can be divided into the "first stimulus" (stimEpoch = 1) and the "second stimulus" (stimEpoch = 2), enabling peak maximum determination within a small window of the trace.
#  The epochs are defined with respect to the stimuli found in stimList. 
#  get_stimEpochs.R accepts the following:
#  df, a dataframe
#  keysEpoch, the collection of column variables that will be combined to define the stimEpoch for each ROI_key. This is potentially deprecated.
#  groupers, the collection of column variables that are used to distinguish each trace from one another. As stimEpochs are mutated into the dataframe, this does not eliminate previously created columns from mutate_condition.





					### adding stimEpochs to dataframe
					df = df_mutated
					rm(df_mutated)
					keysEpoch = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber",'stimEpoch')
					groupers = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol", "exposeNum", "ROINumber") 
					source(paste0(path,"peakAnalysis/get_stimEpochs.R") )

					



					### adding trackROI_key to dataframe. This is necessary to track bouton function across experiments (replicates, fluid exchange rounds, stimulus protocols)
					df=df_stimEpoch
					rm(df_stimEpoch)
					some_conditions = c("sensor","dish","plate", "region","ROINumber" )
					other_conditions=NULL
					.varname1 = "trackROI_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1, df_stimEpoch)




					### adding vid_key to dataframe. This allows measurement of the number of unique recordings for each dataset.
					df = df_mutated
					rm(df_mutated)
					some_conditions = c("dish","plate", "region")
					other_conditions=NULL
					.varname1 = "vid_key" 
					source(paste0(path,"peakAnalysis/get_mutate_condition.R"))
					rm(df,some_conditions,other_conditions,.varname1)


## expFitter.R is a sub-routine that measures the exponential decay function of each peak in the dataset. 
#  expFitter.R relies on nls_multstart, a nonlinear least squares fitting routine that allows arbitrary definition of the function to be fit.
#  Here, we use nls_multstart to perform a group-wise, single-exponential decay fit on thousands of individual dF/F vs. time series.
#  peakIDs, which define a span of points around identified iGluSnFR3 transients, along with other column groupers, are used to distinguish peaks from one another.
#  Starting from the maximum point of each peak and continuing until the peakID terminates in time, the single exponential decay function is found for each peakID. 
#  Once the fitting routine is finished, additional parameters are added to the data.frame to facilitate ease of plotting transients using tracePlotter.R

#  Info from the fitting routine, including:
#  the pre-exponential factor, A
#  the decay time constant, tau_decay
#  and the residuals
#  are stored in a nested data.frame format. Subsequent scripts unpack the nested data.frame to extract tau_decay and A for each peak in the dataset. 

####### for reference, the dataset used in Figure 4 requires ~30 minutes to generate all of the exponential decay fits. #######




					### exponential fitting across the dataset.
					df = df_mutated 
					rm(df_mutated)
					groupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","exposeNum",  "ROINumber",'peakID') 	
					source(paste0(path,"peakAnalysis/expFitter.R"))
					df = fitPeaks
					rm(fitPeaks)





## get_peakStats_activity_vs_marker.R is a sub-routine that is called in multiple annotation packages.
#  get_peakstats_activity_vs_marker.R calculates peak statistics for individual peaks in the dataset, including:
#  amplitude, the maximum dF/F for individual peaks
#  tau_decay_ms, the decay time constant in milliseconds of each peak's exponential decay fit
#  interSpike_ms, the time interval between the peak maximum and the preceding electrical stimulus (if applicable)
#  t_rise, the interpolated 10-90% rise time 
#  t_decay, the interpolated 90-10% decay time
#  t_half, the interpolated full-width at half-maximum







					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key","trackROI_key","ROINumber",'peakID')
					subgroupers = c("sensor","marker","segmentation","dish","plate", "region", "Ca", "protocol","vid_key","trackROI_key") 		 
					drops = c('data', '.resid','fit',"(weights)",".fitted")
					keys = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum" )
					tic()
					source(paste0(path,"peakAnalysis/get_peakStats_activity_vs_marker.R"))
          			rm(drops)
  					rm(df_mutated)






### uncomment this section to plot traces for visualization
## tracePlotter.R generates publication quality visualizations of individual traces. 
#  tracePlotter.R accepts:
#  ROI_IDs, the list of column headers to differentiate between traces. Could be deprecated, as ROI_key should work just as well. 
#  n, the factor by which to divide the total number of ROIs. n=100 would result in 1% of the dataset being saved out as .jpeg and .emf at 600 dpi. 
#  hide_titles, a Boolean which determines whether the ROI_key of origin will be displayed as part of the trace.
#  add_max,  a Boolean which adds the peak maximum.
#  add_time_dots, a Boolean which adds all of the additional dots to the peak to define t_half, t_rise, t_decay



				    # ROI_IDs = c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol",  "ROINumber")
					# n=100 #factor by which we should divide the ROIs_sample - e.g. do you want to see  1/5 or 1/10 of the dataset?  
					# pub_quality=TRUE
					# hide_titles=TRUE
					# add_max = TRUE
					# add_time_dots = TRUE
					# source(paste0(path,"peakAnalysis/tracePlotter.R")) 
					# rm(ROI_IDs,add_max,hide_titles,n)

						









  				source(paste0(path,"peakAnalysis/read_multiple_csv.R"))
  					spont_data_dir = "Y:\\Sam/paper1_datasets/spont_v1/Figures/spont_figs_v3_SI/csv_output"
  					setwd(spont_data_dir)
  					prefix = "clean_stats."
  					fname = "spont_phys_v1_clean_stats.csv"

  					spont_stats<- csv_reader(fname = fname, prefix = prefix)

					  prefix = "clean_df."
  					fname = "spont_phys_v1_clean_df.csv"
  					spont_df<- csv_reader(fname = fname, prefix = prefix)






  					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					
					
					source(paste0(path,"peakAnalysis/get_plot_peakStats_spont_evoked_v1.R"))
					
					file_prefix = "spont_vs_evoked_phys_v2_"
					source(paste0(path,"peakAnalysis/get_written_output_csv.R"))
					
					spont_evoked_output<-output_data

													#file_prefix = "spont_phys_v1_"

					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					spont_avg_df = spont_evoked_output$Averages %>% dplyr::filter(protocol_fixed == "Spontaneous")
			    
					ROIs_to_save = c("GluSnFR3-dish05-plate09-region02-ROI0017", #ROI0016
					          "GluSnFR3-dish05-plate09-region02-ROI0021") #ROI00

					source(paste0(path, "peakAnalysis/get_MV_analysis.R"))
					file_prefix = "MV_analysis_"
					output_data = MV_output
					source(paste0(path,"peakAnalysis/get_written_output_csv.R"))


                          # placeholder = NULL
						# 							#empty_plot_area = NULL
						# 			 				margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
						# 							gs = list(avg_peaks, 	#1
						# 							             nostats_scatter_amplitude,
						# 							             histo_amplitude,
						# 							             histo_tau_decay_ms,
						# 							             #histo_t_half,
						# 							             quanta_conv,		#3
						# 							             N_sites_panelF,
						# 							             N_sites,
						# 							             RRP_proportion_distribution

						# 							             )	#6
													

						# 							 hlay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
						# 							               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
						# 							               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
						# 							               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), #5x15
						# 							               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
						# 							               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),


						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
						# 							               c(2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4), #7x16
													               
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8), #7x16
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8), #7x16
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8), #7x16
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8), #7x16
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8), #7x16
						# 							               c(5,5,5,5,5,6,6,6,6,6,7,7,8,8,8,8,8)
						# 							               )
												


												
						# 							 scalefactor=0.7#0.75
								    				
						# 							Figure3 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
						# 						        ggsave(filename="Figure3_v6_Histos_spont_vs_evoked_O-phys_paper.emf",plot=Figure3,dpi=900, units="in",width=34*scalefactor,height=36*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
						# 						        ggsave(filename="Figure3_v6_Histos_spont_vs_evoked_O-phys_paper.jpeg",plot=Figure3,dpi=600, units="in",width=34*scalefactor,height=36*scalefactor, device = "jpeg")
                        							
												     margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))              
													gs = list(avg_peaks, 	#1
													             nostats_scatter_amplitude,
													             quanta_conv,		#3
						 							             N_sites_panelF,
						 							             N_sites,
						 							             RRP_proportion_distribution

													             		#3
													             )	#6
													

													 hlay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
													               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
													               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
													               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2), #5x15
													               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2),
													               c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2),

																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6),
																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6),
																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6),
																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6),
																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6),
																   c(3,3,3,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6,6)


													               )






  													Panel2_Figure1_SFN = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												        ggsave(filename="Panel2_Figure1_SFN.emf",plot=Panel2_Figure1_SFN,dpi=900, units="in",width=40*scalefactor,height=24*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												        ggsave(filename="Panel2_Figure1_SFN.jpeg",plot=Panel2_Figure1_SFN,dpi=900, units="in",width=40*scalefactor,height=24*scalefactor, device = "jpeg")
                        
  
												        
												        
						





  					
			
          #df_traces = df
          #clean_stats = output_data$clean_stats %>% dplyr::filter(protocol_fixed == "Evoked")
				  spont_avg_df = spont_evoked_output$Averages %>% dplyr::filter(protocol_fixed == "Spontaneous")
					groupers =c("sensor","marker","segmentation","TTL_start","dish","plate", "region", "Ca", "protocol","exposeNum","vid_key","ROI_key",  "ROINumber", 'windowedPeakID')
					plotBy = c('Ca')
		 			levels = c('0pt5Ca','1Ca','2Ca','4Ca')
					color_override = c('#7AD151FF','#22A884FF','#414487FF','#440154FF')
					drops = c('data', '.resid','fit','tidied',"(weights)",".fitted")
					ROIs_to_save = c(#"GluSnFR3-dish05-plate09-region02-ROI0021",
										"GluSnFR3-dish05-plate09-region02-ROI0021",
										"GluSnFR3-dish05-plate09-region02-ROI0023")#, #ROI0016
					          #"GluSnFR3-dish05-plate09-region02-ROI0021") #ROI0021
					          #"GluSnFR3-dish05-plate04-region02-ROI0033")#,

					          #c("GluSnFR3-dish05-plate03-region01-ROI0011",
								#	"GluSnFR3-dish05-plate03-region01-ROI0021")

					          #"GluSnFR3-dish05-plate09-region02-ROI0023",
										#"GluSnFR3-dish05-plate09-region02-ROI0024")#,
										# "GluSnFR3-dish05-plate09-region02-ROI0001",
										# "GluSnFR3-dish05-plate09-region02-ROI0002",
										# "GluSnFR3-dish05-plate09-region02-ROI0003",
										# "GluSnFR3-dish05-plate09-region02-ROI0004",
										# "GluSnFR3-dish05-plate09-region02-ROI0005",
										# "GluSnFR3-dish05-plate09-region02-ROI0006",
										# "GluSnFR3-dish05-plate09-region02-ROI0007",
										# "GluSnFR3-dish05-plate09-region02-ROI0008",
										# "GluSnFR3-dish05-plate09-region02-ROI0009",
										# "GluSnFR3-dish05-plate09-region02-ROI0010",
										# "GluSnFR3-dish05-plate09-region02-ROI0011",
										# "GluSnFR3-dish05-plate09-region02-ROI0012",
										# "GluSnFR3-dish05-plate09-region02-ROI0013",
										# "GluSnFR3-dish05-plate09-region02-ROI0014",
										# "GluSnFR3-dish05-plate09-region02-ROI0015",
										# "GluSnFR3-dish05-plate09-region02-ROI0016",
										# "GluSnFR3-dish05-plate09-region02-ROI0017",
										# "GluSnFR3-dish05-plate09-region02-ROI0019",
										# "GluSnFR3-dish05-plate09-region02-ROI0020",
										# "GluSnFR3-dish05-plate09-region02-ROI0022",
										# "GluSnFR3-dish05-plate09-region02-ROI0023",
										# "GluSnFR3-dish05-plate09-region02-ROI0024",
										# "GluSnFR3-dish05-plate09-region02-ROI0025",
										# "GluSnFR3-dish05-plate09-region02-ROI0026",
										# "GluSnFR3-dish05-plate09-region02-ROI0027",
										# "GluSnFR3-dish05-plate09-region02-ROI0028",
										# "GluSnFR3-dish05-plate09-region02-ROI0029",
										# "GluSnFR3-dish05-plate09-region02-ROI0030",
										# "GluSnFR3-dish05-plate09-region02-ROI0031",
										# "GluSnFR3-dish05-plate09-region02-ROI0032",
										# "GluSnFR3-dish05-plate09-region02-ROI0033",
										# "GluSnFR3-dish05-plate09-region02-ROI0034",
										# "GluSnFR3-dish05-plate09-region02-ROI0035") #44 #22
				            #c("GluSnFR3-dish05-plate01-region02-ROI0013","GluSnFR3-dish05-plate01-region02-ROI0017")
					source(paste0(path,"peakAnalysis/get_releaseProb_v2.R"))
					
			

					
					placeholder = NULL
													
									 				margin = theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))              
													
													gs = list(map_output_1__clusterID,
													          
																map_output_4__mean_quanta,
													          
													          tmp_tracePlot1,
													          tmp_tracePlot2,
													          
													          tmp_releaseProb1,
													          tmp_amplitude1,
													          #tmp_N_sites1,

															  tmp_releaseProb2,
													          tmp_amplitude2
													          #tmp_N_sites2

													           )	#14
													

													 hlay <- rbind(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),
													 			   

													 			   c(3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4),
													 			   c(3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4),
													 			   c(3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4),
													 			   c(3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4),
													 			   c(3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4),
													 			   

													 			   c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
													 			   c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
													 			   c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8),
													 			   c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8)
													 			   
																	)
												
												
													 scalefactor=0.75
								    				
													Figure4 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												        ggsave(filename="SI_FIgure2_case_study_.emf",plot=Figure4,dpi=900, units="in",width=32*scalefactor,height=34*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												        ggsave(filename="SI_FIgure2_case_study_.jpeg",plot=Figure4,dpi=600, units="in",width=32*scalefactor,height=34*scalefactor, device = "jpeg")
												  

					placeholder = NULL
													
									 				margin = theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))              
													
													gs = list(peak_probability, 	#1
													             CV_amplitude,	#3
													             CV_tau_decay,	#4
													             CV_dt,
													             
													          clustered_UMAP,
													             clustered_PiGlu,
													             clustered_quanta,   
													             clustered_N_sites,
													             
													          #placeholder,
													          map_output_1__clusterID,
													          
													          tmp_tracePlot1,
													          tmp_releaseProb1,
													          tmp_amplitude1,
													          tmp_N_sites1
													           )	#14
													

													 hlay <- rbind(c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4),
													 			   c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4),
													 			   c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4),
													 			   c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4), #4x14
													 			   
													               
													               c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8),
													               c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8),
													               c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8),
													               c(5,5,5,5,6,6,6,6,7,7,7,7,8,8,8),
													 			         
													 			   c(9,9,9,9,9,9,10,10,10,10,10,10,10,10,10),
													 			   c(9,9,9,9,9,9,10,10,10,10,10,10,10,10,10),	
													 			   c(9,9,9,9,9,9,10,10,10,10,10,10,10,10,10),	
													 			   c(9,9,9,9,9,9,10,10,10,10,10,10,10,10,10),	
													 			   c(9,9,9,9,9,9,11,11,11,12,12,12,13,13,13),
													 			   c(9,9,9,9,9,9,11,11,11,12,12,12,13,13,13),													 			         
													 			   c(NA,NA,NA,NA,NA,NA,11,11,11,12,12,12,13,13,13),
													 			   c(NA,NA,NA,NA,NA,NA,11,11,11,12,12,12,13,13,13)
													 			   
																	)
												
												
													 scalefactor=0.75
								    				
													Figure4 = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												        ggsave(filename="Figure4_v8_physiology_from_tracking_O-phys_paper.emf",plot=Figure4,dpi=900, units="in",width=30*scalefactor,height=32*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												        ggsave(filename="Figure4_v8_physiology_from_tracking_O-phys_paper.jpeg",plot=Figure4,dpi=600, units="in",width=30*scalefactor,height=32*scalefactor, device = "jpeg")
												  


												    margin = theme(plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"))              
													
													gs = list(clustered_avg_traces, 	#1
													             clustered_CV_amplitude,	#3
													             clustered_CV_tau_decay,	#4
													             clustered_CV_dt,
													             
													          clustered_tau_decay,
													             clustered_dt,
													          clustered_t_half,   
													          
													          clustered_t_rise,
													             clustered_t_decay,
													             clustered_Pr,
													             clustered_RRP_ratio,
													             clustered_Q
													             
													           )	#14
													

													 hlay <- rbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
													 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
													 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
													 			   c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), #4x14
													 			   
													               
													               c(NA,NA,NA,2,2,2,3,3,3,4,4,4,NA,NA,NA),
													               c(NA,NA,NA,2,2,2,3,3,3,4,4,4,NA,NA,NA),
													               c(NA,NA,NA,2,2,2,3,3,3,4,4,4,NA,NA,NA),

													               c(5,5,5,6,6,6,7,7,7,8,8,8,9,9,9),
													               c(5,5,5,6,6,6,7,7,7,8,8,8,9,9,9),
													               c(5,5,5,6,6,6,7,7,7,8,8,8,9,9,9),

													               c(NA,10,10,10,10,11,11,11,11,12,12,12,12,NA,NA),
													               c(NA,10,10,10,10,11,11,11,11,12,12,12,12,NA,NA),
													               c(NA,10,10,10,10,11,11,11,11,12,12,12,12,NA,NA),
													               c(NA,10,10,10,10,11,11,11,11,12,12,12,12,NA,NA)   
													 			   
																	)






												     SI_cluster_data = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												        ggsave(filename="SI_clustered_data_extended_fig4.emf",plot=SI_cluster_data,dpi=900, units="in",width=30*scalefactor,height=26*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												        ggsave(filename="SI_clustered_data_extended_fig4.jpeg",plot=SI_cluster_data,dpi=600, units="in",width=30*scalefactor,height=26*scalefactor, device = "jpeg")
												  




												  # file_prefix = "releaseProb_analysis_"
												  # source(paste0(path,"peakAnalysis/get_written_output_csv.R"))
												        

												  #   margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
													
												# 	gs = list( clustered_UMAP,
												# 	             clustered_PiGlu,
												# 	             clustered_quanta,   
												# 	             clustered_N_sites,
												# 	          map_output_1__clusterID
													
												# 	             )	#14
													

												# 	 hlay <- rbind(c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5),
												# 	               c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5),
												# 	               c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5),
												# 	               c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5),
												# 	               c(1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5)

												# 	               									 #8x14 remains
												# 	               # c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	               # c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	               # c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	               # c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
													               
												# 	               # c(7,7,7,7,7,7,7,8,8,8,8,9,9,9,9),
												# 	               # c(7,7,7,7,7,7,7,8,8,8,8,9,9,9,9),
												# 	               # c(7,7,7,7,7,7,7,8,8,8,8,9,9,9,9),
												# 	               # c(7,7,7,7,7,7,7,8,8,8,8,9,9,9,9)
													               
												# 	               )
												
												
												# 	 scalefactor=0.75
								    				
												# 	Panel2_Figure3_partA_SFN = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												  #       ggsave(filename="Panel2_Figure3_partA_SFN.emf",plot=Panel2_Figure3_partA_SFN,dpi=900, units="in",width=44*scalefactor,height=10*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												  #       ggsave(filename="Panel2_Figure3_partA_SFN.jpeg",plot=Panel2_Figure3_partA_SFN,dpi=1200, units="in",width=44*scalefactor,height=10*scalefactor, device = "jpeg")



												  #   # gs = list(#map_output_1__CV_amplitude,
												# 	#            #  map_output_4__mean_quanta,
												# 	#             # map_output_0.5__releaseProbability_perROI#,#8
												# 	#              #map_output_4__release_ratio,

												# 	#               tmp_tracePlot1, #9
												# 	#               tmp_releaseProb1,#10
												# 	#               tmp_amplitude1,#11
													             
												# 	#               tmp_tracePlot2, #12
												# 	#               tmp_releaseProb2,#13
												# 	#               tmp_amplitude2
												# 	#              )	#14
													

												# 	#  hlay <- rbind(#c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3),
												# 	#                #c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3),
												# 	#                #c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3),
												# 	#                #c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3),
												# 	#                #c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3)#,

												# 	#                									 #8x14 remains
												# 	#                c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3),
												# 	#                c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3),
												# 	#                c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3),
												# 	#                c(1,1,1,1,1,1,1,2,2,2,2,3,3,3,3),
													               
												# 	#                c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	#                c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	#                c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6),
												# 	#                c(4,4,4,4,4,4,4,5,5,5,5,6,6,6,6)
													               
												# 	#                )
												
												
												# 	#  scalefactor=0.75
								    				
												# 	# Panel2_Figure2_partB_BRAIN = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
												  #   #     ggsave(filename="Panel2_Figure2_partB_BRAIN.emf",plot=Panel2_Figure2_partB_BRAIN,dpi=900, units="in",width=30*scalefactor,height=16*scalefactor, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
												  #   #     ggsave(filename="Panel2_Figure2_partB_BRAIN.jpeg",plot=Panel2_Figure2_partB_BRAIN,dpi=1200, units="in",width=30*scalefactor,height=16*scalefactor, device = "jpeg")







							file_prefix = "evoked_physiology_v1"
							source(paste0(path,"peakAnalysis/get_written_output_csv.R"))


					
toc()

