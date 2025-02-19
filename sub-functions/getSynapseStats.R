#getSynapseStats.R


tic()
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/interSpike.R") 


cleanInterSpike<-interSpike(df = fitsCleaned, 
								stimList = stimList, 
								groupers = groupers,
								groupers_byPeak = groupers_byPeak, 
								plotBy = plotBy, 
								levels = levels, 
								secondAxis = secondAxis, 
								color_override=color_override, 
								keys = keys,
								peak_keys = peak_keys,
								bin_override = bin_override)


print("Finished getting synapse stats.")
toc()


