
#Written Samuel T Barlow 2.21.2025
#fileManipulate.R is a function wrapper which can coerce Frame-Intensity files from ImageJ into the preferred format for peakFinder.R.
#The script performs two functions:
# 1. the script reads in all the .csv files containing Frame and Intensity data and staples a few expected columns to the data frame. These are: 
#    fileID - the title of the subdirectory from which the file was read
#	 timeStamp - a deprecated column which is likely not needed for full functionality. In the original workflow, timeStamp is the date-time format for camera frame times. 
#	 index - the Frame number, but starting from index 1. 
# 	 dt - an interpolated frame time based on the known exposure time of the recording. This is a simple linear interpolation in the current format (e.g. frame 1 = 0 s, frame 2 = 0.005 s, etc.)

#This script requires the user to specify a few details.
# path - the path to the code repository. please specify the outer shell of the directory (e.g. from Github, Drive:\YourPath\HighThroughputOphys-main)
# NOTE: There is a character limit of ~244 characters in a file path. Paths longer than 244 characters will throw and error as R cannot read the files. Consider renaming your directory to limit the number of characters. 
# dir - the path to your dataset to be processed. This will be the .csv file directory containing the folder "ROIoutputs". 
# NOTE: The script is expecting to find a directory called "ROIoutputs" so that it can recurse into it for all of the .csv files. Do not specify your dir as a subdirectory below "ROIoutputs". 
# frameTime - this should be given in seconds.

#The remainder of the global variables need not be manipulated unless the user has adjusted prior steps in the pipeline.


#uncomment lines 22 and 23 if you have not run this script before
#list_of_packages(c('plyr','dplyr','gridExtra','zoo','tidyr','reshape2','stringr','readr','purrr','tools','tictoc','sjmisc'))
#install.packages(list_of_packages)

#Read in usual dependencies
library(plyr) ###File manipulation
library(dplyr) ###File manipulation
library(gridExtra) ###Tabulating data
library(zoo)  ###not used in this, but mathematical functions
library(tidyr)
library(reshape2)
library(stringr)
library(readr)
library(purrr)
library(tools)
library(tictoc)
library(sjmisc)


#### USER ADJUSTED PARAMETERS ####

path = "Z:\\R Scripts/Github_repo"
dir = "Y:\\Sam/test_pipeline/your_data_in_csv_format" #data directory
frameTime = 0.1#time in seconds that each frame is expected to take. This value is used to interpolate the time points for each indexed frame. 





### PARAMETERS TO TARGET FILES - DO NOT ADJUST UNLESS YOU HAVE MODIFIED OUTPUTS FROM PRIOR IMAGEJ MACROS ###

subfolders = list.dirs(path=dir,full.names=FALSE, recursive=TRUE) 
subfolders = grep("_Results", subfolders, value=TRUE)
pattern1 = "_driftcorr-\\d_Results|_driftcorr_Results|-\\d_Results|_Results"
subfolders_fixed = str_remove(subfolders, pattern = pattern1)

newFilePat = "ROI\\d\\d\\d\\d"
          
tic()
source(paste0(path, "/fileManipulateR/fileRename.R")) ##only run this if you need to rename files still; converts all ROIs to uniform numbering
toc()
