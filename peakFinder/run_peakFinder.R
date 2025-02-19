##run_peakFinder.R


# written by Samuel T Barlow
# 10.3.22
##peakFinder.R is a function-wrapper that will accept a working directory in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 


#changelog 
#12.10.24 - changed window duration from 0.75 to 0.2. toggle_long_peaks = TRUE is best for HNK iGluSnFR3 dataset


path = "Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/"


library(tictoc)

# set working directory
setwd("Y:\\Sam/post_paper1_datasets/HNK/HNK_v2/.csv/HNK_stamped")

segmentedBy = "marker_long_peaks"
TTL_start = 85
showGraphs = FALSE
sensor_override = TRUE
sensor = "GluSnFR3"
groupers =c("sensor","dish","plate", "region","Ca",'protocol',"HNK_phase",'exposeNum',"ROINumber") #"chemical_condition", 
						
#if savePlots is TRUE, need to set a directory and number of plots to save out
savePlots = TRUE
mainDir <- "Y:\\Sam/post_paper1_datasets/HNK/HNK_v2/Figures/peakFinder_outputs_v1" #set directory for Figures
n =100#a number: # of traces / n is the % of plots to save out
# 
tic()
source(paste0(path,"peakFinder/peakFinder.R"))
toc()				