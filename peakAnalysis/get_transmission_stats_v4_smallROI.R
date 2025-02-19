
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/transmission_stats_v4_smallROI.R") 


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
