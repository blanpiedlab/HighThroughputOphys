#get_releaseProb.R



tic()
source(paste0(path,"peakAnalysis/def_stimEpochs.R"))




df_stimEpoch<-def_stimEpochs(df = df, 
								stimList = stimList, 
								groupers=groupers,
								keysEpoch=keysEpoch
								)

print("Finished defining stimulus epochs.")
toc()

