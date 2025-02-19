### get_peakStats_activity_vs_marker.R
### 12.6.23



tic()
source(paste0(path,"peakAnalysis/peakStats_activity_vs_marker.R") )

			
peak_stats<-peakStats(df=df,groupers=groupers,subgroup=subgroupers,list_stimuli=stimList,keys=keys)


print("Finished getting peak stats.")
toc()


