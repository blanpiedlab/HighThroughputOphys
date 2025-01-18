### get_peakStats_activity_vs_marker.R
### 12.6.23



tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/peakStats_spont_phys_v1.R") 

			
peak_stats<-peakStats(df=df,groupers=groupers, remove_cols=drops)


print("Finished getting peak stats.")
toc()


