
source(paste0(path,"peakAnalysis/transmission_stats_v3.R") )


transmissionSummary <- transmission_stats(traces_df = df_traces, 
											interSpike_thresh = NULL, 
											peak_groupers = peak_groupers,
											groupers = groupers,
											positive_ROI_groupers = positive_ROI_groupers, 
											plotBy = plotBy, 
											levels = levels, 
											secondAxis = NULL, 
											color_override=color_override,
											comparisons = my_comparisons)
