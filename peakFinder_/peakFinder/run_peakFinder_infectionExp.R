##run_peakFinder_infectionExp.R


# written by Samuel T Barlow
# 10.10.22
##peakFinder.R is a function-wrapper that will accept a working directory in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 




# set working directory
setwd("Y:\\Sam/2022-09-16/datasets/full_PP")

segmentedBy = "activity"
TTL_start = 199

# 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/peakFinder_infectionExp.R")
						