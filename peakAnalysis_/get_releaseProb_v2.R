
tic()
source(paste0(path, "peakAnalysis/releaseProb_v2.R")) 

#break_me""
			
output_data<-releaseProb_stats(peak_stat_df = peak_stats, df = df, spont_avgs = spont_avg_df, MV_analysis = MV_output, groupers = groupers, color_override = color_override, remove_cols=drops, plotBy = plotBy, levels=levels,remove_col = drops,ROIs_to_save = ROIs_to_save)

print("Finished analyzing synapse heterogeneity.")

toc()