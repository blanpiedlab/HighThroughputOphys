#get_releaseProb.R



tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v4 - markerSegmentation/def_stimEpochs.R") 




sync_async_fitsLabeled_stimKeyed<-def_stimEpochs(df = sync_async_fitsLabeled, 
								stimList = stimList, 
								groupers=groupers,
								keys = keys,
								keysROI=keysROI,
								keysEpoch=keysEpoch
								)


print("Finished defining stimulus epochs.")
toc()

