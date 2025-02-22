##run_peakFinder.R
# written by Samuel T Barlow
# updated 2.20.25

##peakFinder.R is a function-wrapper that will accept a directory of .csv files containing Time and Intensity information in which to identify peaks and output a .csv containing the peak-annotated data frame. 
# the goal of this function-wrapper is to minimize my interactions with the code as much as possible and make data analysis, figure generation easier by making peak identification a single-use event for a given dataset. 

### if you are running this script for the first time, uncomment line 10 and 11. 
#list_of_packages = c('plyr','dplyr','ggplot2','zoo','ggthemes','tidyr','reshape2','stringr','readr','purrr','extrafont','tictoc','broom','gridExtra','gsignal','ggeasy','scales')
#install.packages(list_of_packages)

library(tictoc)

path = "YourDrive:\\HighThroughputOphys-main/"



# set working directory
# Please find the data at the following link: https://drive.google.com/drive/folders/1_Mm8QNhFO6zKyGGvtbehgSkhrEo5kW_H?usp=drive_link
setwd("YourDrive:\\HighThroughputOphys-main/demos/run_peakFinder/demo_data_peakFinder")

segmentedBy = "marker"			# in Barlow et al. bioxRxiv 2024, we segmented videos on the basis of a fluorescent marker (Synaptophysin-mRuby3 expression), or iGluSnFR3 activity itself. These were referred to as "marker" and "activity", respectively. 
TTL_start = 195 			# in Barlow et al. bioRxiv 2024, we use a TTL pulse at a defined frame to initiate stimulus protocols. The frame number can be used to interpolate the true stimulus position in later analyses. 
toggle_short_peaks = FALSE 		# Toggle to TRUE if the fluorescence transients in your data are frequently convoluted with one another due to the frequency of the process. toggle_short_peaks will make it more likely to find multiple local maxima in highly convoluted portions of the trace, but may impact baseline identification negatively.  
showGraphs = FALSE 			# Toggle to TRUE if you'd like to see what the algorithm is up to in real time. Note that this will slow down the computation significantly.  
sensor_override = TRUE 			# Allows the user to set the fluorescent sensor for the dataset. This determines which baseline adjustment method and which signal kinetics will be used the dataset.
sensor = "GluSnFR3" 			# If sensor_override is TRUE, here you can define a sensor from among the following: c("GluSnFR3", "JF646"). JF646 refers to JF646-BAPTA; refer to Deo et al. JACS 2019. Future iterations of the code will be able to select from c('GluSnFR3',"JF646",'GCaMP8m','GCaMP8f','GCaMP6f','jRGECO').  At present, if you wish to perform calcium imaging analysis, select 'JF646'.

groupers =c("sensor","dish","plate", "region","Ca",'protocol','exposeNum',"ROINumber") # a list of column names with string identifiers used to ensure only one trace is under consideration at a time. Implementation with basic identifiers coming soon.   


### if savePlots is TRUE, need to set a directory and number of plots to save out
savePlots = TRUE
mainDir <- "YourDrive:\\HighThroughputOphys-main/demos/run_peakFinder/Figures_dir" #set directory for Figures
n =10 # a number. Analogous to a dilution factor, this is the (total number of traces) / n. So n=1 would give you 100% of the traces as plots, while n=100 would give you 1% of the traces as plots. 



tic()
source(paste0(path,"peakFinder/peakFinder.R"))
toc()				
