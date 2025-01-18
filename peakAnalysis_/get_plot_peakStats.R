
tic()
source(paste0(path, "peakAnalysis/plot_peakStats_spont_activity_vs_marker.R") ) 

#break_me""
			
output_data<-plot_peakStats(peak_stat_df = peak_stats, df = df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col=drops)

print("Finished getting peak stats by segmentation.")

toc()

