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


#if dirname sort is set to TRUE, the directory of the zstack should be used to regex match for ROImaps



####Set the directory manually for extracting data####
workingDir = "Y:\\Sam/post_paper1_datasets/HNK/initial_exp_1Ca_30uM_HNK/imageFiles" # this is probably redundant
setwd(workingDir)
stack_dir = "Y:\\Sam/post_paper1_datasets/HNK/initial_exp_1Ca_30uM_HNK/imageFiles/zstackOutputs_sliced_SQ"

vid_dir = "Y:\\Sam/post_paper1_datasets/HNK/initial_exp_1Ca_30uM_HNK/imageFiles/toBeAnalyzed/drift_corrected_tiffs"

stack_pat = "ROImap.zip"#"C1(.*?).tif"#  #

#match_pats is deprecated. Now using vid_regex matching between GluSnFR3 and the Calcium identity which includes all plate-region identifiers. This seems to work well
#match_pats = "_Results_\\d|_Results_\\d\\d" 


vid_regex = "dish(.*?)region\\d\\d_"
  
  
  ##"^\\(GluSnFR3_)(PSD95FingR|SyPhy_mRuby3)(_dendrite|_axon)(_05mMCa|_1mMCa|_4mMCa)(_dish\\d)(_plate\\d)(_region\\d).*\\.tif$" 



source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/fileManipulateR/copy_ROImaps_v2.R")


tic()

copy_ROImaps(zstack_dir = stack_dir, zstack_pat = stack_pat, video_dir = vid_dir, video_pat=vid_regex, matches = vid_regex, dirname_sort = TRUE)

toc()