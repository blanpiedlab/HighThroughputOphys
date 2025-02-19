
tic()
source(paste0(path, "peakAnalysis/MV_analysis.R")) 

#break_me""
			
MV_output<-MV_analysis(peak_stat_df = peak_stats, df = df, spont_avgs = spont_avg_df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops,ROIs_to_save = ROIs_to_save)

print("Finished getting peak stats by segmentation.")

toc()