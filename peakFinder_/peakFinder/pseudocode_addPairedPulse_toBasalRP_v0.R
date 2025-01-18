##extract_BasalRP_from_PairedPulse.R


# written by Samuel T Barlow
# 10.27.22
# extract_BasalRP_from_PairedPulse.R is an optional script designed to executed when peakFinder is applied to a Paired Pulse dataset. If extract_BasalRP = TRUE, extract_BasalRP_from_PairedPulse.R will execute as an auxiliary process following peakFinder_Lite.R

# extract_BasalRP_from_PairedPulse.R accepts two inputs: 
# 1. a dataframe containing our PairedPulse data
# 2. a vector containing the protocols to extract - ideally these will also be presented in their preferred sequencing , e.g. PP1000 x3, PP750 x3, PP500 x3, PP250 x3... 

# extract BasalRP_from_PairedPulse.R will then produce a single output - a modified version of the Paired Pulse dataframe which can be read into a BasalRP analysis module as a separate annotated_csv. Ideally, the output will be such that the operator can differentiate 
# that the annotated_csv came from a PairedPulse dataset, but it should not be differentiable from the perspective of downstream BasalRP analyses.

#The overall objective of this script is to add replicates of BasalRP to each ROI. 
# Since the protocols are run in a sequence, this means that each ROI should have :
# 5 replicates from BasalRP
# 3 replicates from each Paired pulse set in which the peaks are separable (this dataset, 250,500,750,1000)
# total : 5*BasalRP + 3*(4 Paired Pulse) = 17 replicates.
# This will give us a much finer-grained analysis of release Probability. 



# The approach: First, we will extract only the "protocol" that we want. This will be defined by an input vector containing the desired protocol strings
# df %>% dplyr::filter(protocol %in% protocol_list)

# Once we have the protocol list, it gets a little tricky - we need to have ROIkey 








#after talking to Poorna - consider creating a master-key of trackROIs with their number of replicates from BasalRP. 

# from BasalRP



# read in master key.
# generate trackROIs in paired pulse.
# create columns - total replicates in BasalRP, total_PP1000, total_PP750, total_PP500, totalPP250
# create a dictionary of FROM: protocol_repl labels
# create a dictionary of TO: protocol_repl labels
# left_join original dataframe +trackROIs with TO: protocol_repl labels
# drop existing replicate column
# drop existing protocol labels
# TO: protocol = BasalRP
# TO: exposeNum = new replicate names (paste0('repl',value of true-repl column))
# drop newExposeNum == exposeNum
# drop newExposeNum
# drop trackROIs
# output to annotated.csv






### set annotated BasalRP directory
setwd()
trackROI_keys = c("sensor","marker","neuronSegment","segmentation","TTL_start","dish","plate", "region", "Ca", "ROINumber")
groupers=c()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/annotated_csv_reader.R")

#read in the BasalRP dataset, get trackROI_keys, define existing replicates per trackedROI
#df %>% group_by_at(groupers)


source("getBasalRP_replicates")
rm(BasalRP_dataset)
### set the PP annotated csv directory

setwd()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/annotated_csv_reader.R")

#read in the PP dataset, get trackROI_keys, define existing replicates per protocol (group_by_at(groupers)summarise())

source("getExisting_PPreplicates")

rm(PP_dataset) #memory


source("join_datasets")
### combine datasets and get 4 sums - BasalRP, BasalRP+protocol1, sum + protocol2, sum +protocol3....



### need to retain protocols for replicate_summer so it can dplyr::group_by_at("protocol","trackROI_keys")
### feed protocol list #### needs more thought to execute
protocol_order = c('PP1000', 'PP750','PP500','PP250')

source("replicate_summer")

### extract replicate number gsub() 


### new_replicate = case_when( current_protocol == protocol_order[1] ~ paste0(BasalRP_replicates+ gsub(exposeNum)),
													#current_protocol == protocol_order[2]....)

#drop extra sums

### re-read Paired Pulse annotated csv, left_join via trackROIs, drop rename old_protocol to protocol == "BasalRP", rename old_exposeNum to exposeNum == paste0("repl",toReplicate), get TTL_time, chop all absoluteTimes past TTL_time + 220 ms, drop TTL_time
setwd()
source("PPtoBasalRP_annotated_csv_output")
## output as a modified title so we know it came from PP datasets



# set a working directory


# set working directory
dataDir = "Y:\\Sam/2022-10-23/Figures/full_PP_v0/annotated_csv"

setwd(dataDir)
mainDir <- "Y:\\Sam/2022-10-23/Figures/full_PP_v0" #set directory for Figures




# 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/peakAnalysis.R")
						