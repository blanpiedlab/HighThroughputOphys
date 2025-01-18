
tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/PPR_plotter_v4_Figure4.R") 

#break_me""
			
PPR_stats<-PPR_stats(df=df,groupers=groupers,stimEpoch_groupers=stimEpoch_groupers,plotBy=plotBy,levels=levels, color_override = color_override,keys = keys,save_ROIs = ROIs_to_save, keep_PP = keep_PP)

print("Finished getting pairedPulse ratios.")

toc()