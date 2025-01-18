
tic()
source(paste0(path,"peakAnalysis/plot_peakStats_spont_evoked_v1.R") )

#break_me""
			
output_data<-plot_peakStats(AP_stat_df = peak_stats, AP_df = df, spont_stat_df = spont_stats, spont_df = spont_df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops)

print("Finished getting peak stats for Figure 3.")

toc()

