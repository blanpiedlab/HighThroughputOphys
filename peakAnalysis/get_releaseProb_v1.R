
tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/releaseProb_v1.R") 

#break_me""
			
output_data<-releaseProb_stats(peak_stat_df = peak_stats, df = df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops,ROIs_to_save = ROIs_to_save)

print("Finished getting peak stats by segmentation.")

toc()