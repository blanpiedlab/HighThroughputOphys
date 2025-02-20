# annotated_csv_output.R

# written by Samuel T Barlow
# 10.3.22

# annotated_csv_output.R accepts a dataframe from peakFinder_Lite.R and creates a .csv of it for future manipulations.
# the objective of this script is to split the analysis pipeline into multiple "master scripts" that are fully generalized and accelerate data analysis timelines.

dir = getwd()
mainDir = dirname(dir)
outerDir = dirname(mainDir)
#as_date_origin = str_extract(outerDir, "\\d\\d\\d\\d-\\d\\d-\\d\\d")
subDir = 'annotated_csv'


if (file.exists(file.path(mainDir,subDir))){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir, subDir))
  
}


#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R")
						



dataset_id = unique(peaks$sensor)
print("Now attempting to output an annotated_csv. Printing the existing dataset_ids.")
print(dataset_id)
print(peaks)

print(paste0("There are ", length(dataset_id), " unique sensor datasets here." ) ) 

for (id in dataset_id ) {

	sensor_subset <- peaks %>% dplyr::filter(sensor == id) 

	protocol_id = unique(sensor_subset$protocol)

	print(paste0("There are ", length(protocol_id), " unique protocol datasets here." ) ) 
	#write.csv(sensor_subset, file = paste(as_date_origin, id, k, ".csv",sep="_"))
	for (k in protocol_id) {

	 	protocol_subset <- sensor_subset %>% dplyr::filter(protocol == k)


	 	#TTL_id = first(unique(protocol_subset$TTL_start) )
	 	segmentation_id = first(unique(protocol_subset$segmentation))
 		#write.csv(protocol_subset, file = paste( id, k,"segmentation-",segmentation_id,"TTL_start_",TTL_id, ".csv",sep="_"))
		write.csv(protocol_subset, file = paste( id, k,"segmentation-",segmentation_id, ".csv",sep="_"))

	 }

}
rm(list=ls())

