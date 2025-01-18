


library(tictoc)
tic()

# set working directory
dataDir = "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/.csv/annotated_csv"

setwd(dataDir)

PP_mainDir <- "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/Figures/PP_v3"
mainDir = PP_mainDir


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/peak_annotation_module_PPedit.R")
toc()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/peak_analysis_module_PairedPulse_v2.R")	