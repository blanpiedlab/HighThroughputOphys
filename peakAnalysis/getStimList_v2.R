#getStimList_v2.R
source(paste0(path,"peakAnalysis/def_stim_paradigms_v4.R")) 

stimList = def_stim_paradigms(dataframe=df, APcount = num_APs,keys=keys,offset=s88x_offset)
