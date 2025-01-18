source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/single_ROI_display.R") 


dual_record_display <- single_ROI_display(ROI_traces = small_df_traces, 
											interSpike_thresh = NULL, 
											peak_groupers = peak_groupers,
											groupers = groupers,
											positive_ROI_groupers = positive_ROI_groupers, 
											plotBy = plotBy, 
											levels = levels, 
											secondAxis = NULL, 
											color_override=color_override,
											comparisons = my_comparisons,
											save_ROIs = ROIs_to_save,
											tag_var = tmp_tag,
											display_transmission = disp_trans)
