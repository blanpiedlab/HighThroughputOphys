#Copy ROImaps to their respective folders

#Read in usual dependencies
library(plyr) ###File manipulation
library(dplyr) ###File manipulation
library(ggplot2) ###plotting
library(gridExtra) ###Tabulating data
library(zoo)  ###not used in this, but mathematical functions
library(ggthemes) ###other ggplot2 graph themes
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



####Set the directory manually for extracting data####
workingDir = "Y:\\Sam/2023-03-07/imageFiles/.ims/512_split_2"
setwd(workingDir)
stack_dir = "Y:\\Sam/2022-12-12/imageFiles/markerSegmentation/zstackOutputs"

vid_dir = "Y:\\Sam/2022-12-12/imageFiles/markerSegmentation/toBeAnalyzed-marker"

stack_pat = "ROImap.zip"
match_pats = "GluSnFR3(.*?)region\\d"


vid_pats = c("GluSnFR3", ".tif")
  
  
  ##"^\\(GluSnFR3_)(PSD95FingR|SyPhy_mRuby3)(_dendrite|_axon)(_05mMCa|_1mMCa|_4mMCa)(_dish\\d)(_plate\\d)(_region\\d).*\\.tif$" 



source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/fileManipulateR/copy_ROImaps.R")


tic()

copy_ROImaps(zstack_dir = stack_dir, zstack_pat = stack_pat, video_dir = vid_dir, video_pat=vid_pats, matches = match_pats)

toc()