##peakFinder.R


# written by Samuel T Barlow
# 10.3.22
##peakFinder.R is a function-wrapper that will accept a working directory in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 




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

tic()

#dir = setwd("Y:\\Sam/2022-09-08/datasets/fullPP")
						###
						####

						groupers = c("sensor","chemical_condition", "dish", "plate", "region", "replicate", "ROINumber")							
						#read data - is ready for primetime with whatever data 
						source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/read_csv_filename_mini.R")
						# groupers for the first set of scripts - this shouldn't change much experiment to experiment
						groupers = c("sensor", "chemical_condition","dish", "plate", "region", "replicate", "ROINumber")		

											
						

						#peak identification
						#showGraphs = TRUE 																			# set to TRUE or FALSE to see the graph outputs from peakFinder_Heavy.R/peakFinder_Lite.R
																																		# for a single video dataset with ~35 ROIs and 76 identified peaks, showGraphs = TRUE takes 73 seconds.
																																		# for a single video dataset with ~35 ROIs and 76 identified peaks, showGraphs = FALSE takes 59 seconds.   
						
						source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/flagOutliers_mini.R") 
						
						source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/peakFinder_Heavy_mini.R") 



						#Keys for ROIs write here,
						#These drive peakFinder showGraphs output
						#To maintain with previous versions, will drop column
						ROI_keys = c("sensor", "chemical_condition","dish", "plate", "region", "replicate", "ROINumber")							
						


						## can write in a variance filter a la "peaks that pass"
						## probably groupers works, but working in ROI_keys might be more elegant
						source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/peakFinder_Lite_mini.R") 
						
						
						source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/annotated_csv_output_mini.R") 
						
						

toc()