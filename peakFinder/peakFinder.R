##peakFinder.R


# written by Samuel T Barlow
# 10.3.22
##peakFinder.R is a function-wrapper that will accept a working directory in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 




#necessary libraries
library(plyr) ###File manipulation
library(dplyr) ###File, dataframe manipulation
library(ggplot2) ###plotting
library(zoo)  ###not used in this, but mathematical functions
library(ggthemes) ###other ggplot2 graph themes
library(tidyr) #data frame manipulation
library(reshape2) #data frame manipulation
library(stringr) # string manipulation
library(readr)
library(ggpubr) # higher level ggplot functions, e.g. statistical testing
library(purrr) # dataframe manipulation and nested dataframes
library(extrafont) # extended font selection for ggplot2
library(tictoc) # time benchmarking
library(broom) # dataframe manipulation and nested dataframes
library(gridExtra) # grid layouts for saving multi-figure plots
library("RColorBrewer") # extended color library for ggplot2

#tictoc allows measurement of code run time. 

tic()

### read_csv_filename.R is a function that recurses into a directory to extract all of the .csv files containing intensity-time traces. 
#   users can interact with read_csv_filename.R to generate the desired grouping variables that distinguish individual traces from one another.

#   The basic structure of read_csv_filename.R is to use the package 'stringr' to identify character-based grouping variables from the "source" column
#   	of the intensity-time .csv files. Source contains the entire filestring of the original video file, which are deliberately named to capture all 
#		variable information from the optical physiology recording. Example information that can be found in Source from the experiments in this manuscript are:
#		Sensor: which sensor was used to measure physiology? e.g. GluSnFR3, JF646
#		dish: which 12-well plate did the coverslip come fom? 
#		plate: which coverslip from the 12-well plate do these images correspond to? 
#		region: which imaging region (numeric label) does this video correspond to?
#		Ca: Which calcium concentration (in mM) was present in this video? 
#		exposeNum: what replicate was this recording from?

#	read_csv_filename.R takes advantage of str_detect() along with if/else control statements to identify potential character strings, and if they are detected, add them to the dataframe using str_extract(). 



						source(paste0(path,"/peakFinder/read_csv_filename.R"))
						
																																		# for a single video dataset with ~35 ROIs and 76 identified peaks, showGraphs = TRUE takes 73 seconds.
						

#flag_outliers.R is the first layer of outlier discrimination that leads to the detection of fluorescence transients. 
#flag_outliers.R is simply a median filter that flags all indices which lie outside a rolling median (0.75 s span) +- 1.5σ (σ is standard deviation of the trace).


						ROI_keys = groupers																												# for a single video dataset with ~35 ROIs and 76 identified peaks, showGraphs = FALSE takes 59 seconds.   
						source(paste0(path,"/peakFinder/flagOutliers.R"))
						

#peakFinder_Heavy.R is the second layer of outlier discrimination that leads to the detection of fluorescence transients. 
#peakFinder_Heavy.R takes the rolling median approximation of the baseline from flagOutliers.R, 
#		excludes all indices that were flagged as outliers (coerce to NA), and replaces the NAs with the last non-NA observation (na.locf()).
#		From this new baseline approximation, a rolling average is calculated, and this is treated as the true baseline for each intensity-time trace.
#		Signals are then generously threshold identified using a Schmitt trigger with a lower threshold of 1.5σ and an upper threshold of 3.5σ (here, σ is the standard deviation of ONLY the putative baseline indices). 
#		In this way, all indices associated with a putative fluorescence signal (above the Schmitt trigger) are flagged. Crucially, these flagged indices are padded on either side of the transient to capture the full extend of the signal
#		beyond the limits of the Schmitt trigger identification. This padding function is derived from the published kinetics of each sensor (2.5*the rise time, 2.5* the decay time). 
#		
#		There are edge cases where this padding function leads to unexpected behavior when transients are convoluted with one another (as in JF646 spontaneous recordings), which can lead to a lower rate of signal detection.
#		If toggle_short_peaks = TRUE, the padding function will instead only pad the signal indices  with 0.5*rise time and 0.5*decay time, which can improve signal detection in highly active traces.
#
#		peakFinder_Heavy.R provides a first approximation of dF/F.
#
#		Importantly, peakFinder_Heavy.R is expecting a dataset with only one sensor modality at a time. Trying to run the script with a dataset that contains sensors with profoundly different kinetics will lead to undesirable behavior.


						
						ROI_keys = groupers
						source(paste0(path,"/peakFinder/peakFinder_Heavy.R"))
						
#peakFinder_Lite.R is the final layer that outputs baseline corrected traces in terms of dF/F with peaks annotated for downstream analysis. 
#		peakFinder_Lite.R accepts trace with signal indices found by peakFinder_Heavy.R and outlier indices found by flagOutliers.R and
#		computes a new dF/F using na.locf() and rolling averaging. Peaks are identified in the same way as before, but with a more stringent Schmitt trigger (1.5σ, 5σ).
#		The output of peakFinder_Lite.R is a fully annotated dataframe with found fluorescence transients and dF/F calculations. 



						ROI_keys = groupers
						source(paste0(path,"/peakFinder/peakFinder_Lite.R"))
						

#annotated_csv_output.R is an output function that writes the annotated dataframe to .csv with the specified directory. These are read into R for downstream analyses. 
						
						source(paste0(path,"/peakFinder/annotated_csv_output.R"))
						
						

toc()