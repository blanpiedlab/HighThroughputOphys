#Basic pipeline for file manipulation

#Read in usual dependencies
library(plyr) ###File manipulation
library(dplyr) ###File manipulation
library(ggplot2) ###plotting
library(gridExtra) ###Tabulating data
library(zoo)  ###not used in this, but mathematical functions
library(ggthemes) ###other ggplot2 graph themes
library(DescTools) ###Calculate Area under curve (AUC)
library(tidyr)
library(reshape2)
library(stringr)
library(readr)
library(ggpubr)
library(purrr)
library(extrafont)
library(ggforce)
library(tools)
library(tictoc)
library(sjmisc)

multichannel = FALSE

####Set the directory manually for extracting data####

dir = "Y:\\Sam/tmp/HNK_unstamped"
subfolders = list.dirs(path=dir,full.names=FALSE, recursive=TRUE)
subfolders = grep("_Results", subfolders, value=TRUE)
pattern1 = "_driftcorr-\\d_Results|_driftcorr_Results|-\\d_Results|_Results"
pattern2 = "ROIoutputs/"


###for drift correction handling
#pattern1 = "_driftcorr_Results|_Results"
#pattern2 = "ROIoutputs/"


 if (multichannel == TRUE) {
   sensorID = unique(str_extract(subfolders,pattern=pattern1) ) 
   if (sensorID == "_driftcorr-1_Results" |sensorID == "-1_Results") {
       sensor = "GluSnFR3"
   } else if (sensorID == "_driftcorr-2_Results" |sensorID == "-2_Results") {
       sensor = "JF646"
   }
 }

subfolders_fixed = str_remove(subfolders, pattern = pattern1)
match = unlist(str_remove(subfolders_fixed, pattern = pattern2))

newFilePat = "ROI\\d\\d\\d\\d"
          
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/fileManipulateR/fileRename.R") ##only run this if you need to rename files still; converts all ROIs to uniform numbering

#removed this on 11.229.21: edit 11.19.21 - subset stampfile[200:1200,] to extract only the frames that were sliced out of spontaneous recordings.
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/fileManipulateR/timeStampR.R") #creates a new directory for all timelapses
