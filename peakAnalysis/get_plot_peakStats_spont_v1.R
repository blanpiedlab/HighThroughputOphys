
tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/plot_peakStats_spont_v1.R") 

#break_me""
			
output_data<-plot_peakStats(peak_stat_df = peak_stats, df = df, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col=drops)

print("Finished getting peak stats from spontaneous data.")

toc()

