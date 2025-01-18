#getReplicates.R
#
# Written by Samuel T Barlow 10/24/22
#
# getReplicates.R applies the function find_replicates to the specified data.frame and outputs a revised data.frame including the column "replicates"
# this will improve the flexibility of downstream metrics that divide by the number of replicates (e.g. releaseProbability)

source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/find_replicates.R")
			
df_mutated<- find_replicates(df = df, groupers=groupers, colname_toFind=colname_toFind) 
