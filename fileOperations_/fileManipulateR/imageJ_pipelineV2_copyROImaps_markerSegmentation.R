#Copy ROImaps to their respective folders

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



####Set the directory manually for extracting data####

stack_dir = "Y:\\Sam/2022-09-08/imageFiles/test_Rcode/zstackOutputs"
stack_pat = "ROImap.zip"

vid_dir = "Y:\\Sam/2022-09-08/imageFiles/test_Rcode/test_findRegex"
vid_pat = "^\\(GluSnFR3_)(PSD95FingR|SyPhy_mRuby3)(_dendrite|_axon)(_05mMCa|_1mMCa|_4mMCa)(_dish\\d)(_plate\\d)(_region\\d).*\\.tif$") 



source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/fileManipulateR/copy_ROImaps.R")


tic()

copy_ROImaps(zstack_dir = stack_dir, zstack_pat = stack_pat, video_dir = vid_dir, video_pat=vid_pat)

toc()