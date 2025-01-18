
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/def_peak_halfwidths_rise_decaytimes_v2_spont.R") 

#break_me""
			
tic()

output_data<-get_peak_halfwidths_rise_decaytimes_v1(df = df, groupers = groupers,remove_cols=drops)

print("Finished getting peak stats by segmentation.")

toc()

