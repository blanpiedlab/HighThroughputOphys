##run_mini_peakFinder.R


# written by Samuel T Barlow
# 10.3.22
##peakFinder.R is a function-wrapper that will accept a working directory in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 


#clear all vars
rm(list=ls())


# set working directory
setwd("Y:\\Sam/APVdataset/.csv/stamped/JF646")

#segmentedBy = "marker"
#TTL_start = 195
showGraphs = FALSE

#if savePlots is TRUE, need to set a directory and number of plots to save out
savePlots = TRUE
mainDir <- "Y:\\Sam/APVdataset/Figures/mini_peakFinder_outputs_JF646_2" #set directory for Figures
n = 1#a number: # of traces / n is the % of plots to save out
# 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakFinder/peakFinder_mini.R")
						