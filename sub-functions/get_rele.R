#get_releaseProb.R

source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/def_releaseProb.R") 


keys = c('plate','region','marker','neuronSegment','stimParadigm','exposeNum','ROINumber')
groupersStimuli = c('neuronSegment',"laserPower", "exposureTime", "plate", "region", "Ca_str", "timeClass", "stimParadigm","marker", "exposeNum", "ROINumber",'stimEpoch',"peakID")

groupersROI = c('neuronSegment',"laserPower", "exposureTime", "plate", "region", "Ca_str", "stimParadigm","marker", "exposeNum", "ROINumber")

releaseStats_perStimulus = def_releaseProb(df=fitsCleaned, groupersStimuli=groupersStimuli, groupersROI=groupersROI,keys=keys)