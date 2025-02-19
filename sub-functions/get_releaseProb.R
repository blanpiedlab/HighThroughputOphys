#get_releaseProb.R
subDir <- "cumulativePlots"

source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 

source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v4 - markerSegmentation/def_releaseProb.R") 



releaseStats_perStimulus = def_releaseProb(df_all=sync_async_fitsLabeled_stimKeyed,df_clean=fitsCleaned, 
														groupers_byPeak=groupers_byPeak, 
														groupers_byEpoch=groupers_byEpoch,
														groupers_byROI = groupers_byROI,
														#keysROI=keysROI,
														#keysEpoch=keysEpoch,
														plotBy = plotBy,
														levels = levels,
														varList = varList,
														secondAxis = secondAxis,
														all_ROIs_override = all_ROIs_override)