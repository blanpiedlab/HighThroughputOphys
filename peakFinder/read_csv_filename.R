
#Read in library dependencies
library(plyr) ###File manipulation
library(dplyr) ###File manipulation
library(ggplot2) ###plotting
library(zoo)  ###not used in this, but mathematical functions
library(ggthemes) ###other ggplot2 graph themes
library(tidyr)
library(reshape2)
library(stringr)
library(readr)
library(ggpubr)
library(purrr)
library(extrafont)
library(ggforce)
library(tictoc)


## Read files into R 

filenames<- list.files(pattern=".csv", full.names=TRUE, recursive=TRUE) #Recurse through all subdirectories in the directory


#Define read_csv_filename, which names existing columns and adds new columns according to extracted strings
read_csv_filename <- function(filename, TTL = TTL_start, segmentationMethod = segmentedBy){
		  ret <- read.csv(filename)
						colnames(ret)[1]<- "frameNumber"
					  colnames(ret)[2]<- "intensity"
						colnames(ret)[3]<- "fileID"
						colnames(ret)[4]<- "index"
						colnames(ret)[5]<- "timeStamp"
						colnames(ret)[6]<- "absoluteTime"
						
						if("sensor" %in% colnames(ret)) {
							#print( "yup, there's already a column #7.")
							colnames(ret)[7]<- "sensor"
						} else if(sensor_override == TRUE) {
							ret$sensor = sensor
						} else {
							ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|JF646|GCaMP8f|GCaMP8m")
						}

						ret$Source <- filename

						if( unique(str_detect(ret$fileID, "nomarkers|SyPhy-mRuby3|PSD95FingR|Syn1a|SyPhy"))) {
							ret$marker<-str_extract(ret$fileID, "nomarkers|SyPhy-mRuby3|PSD95FingR|Syn1a|SyPhy")
						}								
						ret$dish<- str_extract(ret$fileID, "dish\\d\\d|dish\\d")
						ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")
						ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
						

						if(unique(str_detect(ret$fileID, "axon|dendrite") ) ) {
							ret$neuronSegment<- str_extract(ret$fileID, "axon|dendrite")
						}

						if(unique(str_detect(ret$fileID, "spont|singleAP|PP\\d\\d\\d|PP\\d\\d") ) ) {
							ret$protocol<- str_extract(ret$fileID, "spont|singleAP|PP\\d\\d\\d|PP\\d\\d")
						} else {
							ret$protocol<- "spont"	
						}
																				#			#				#			#		#
						if(unique(str_detect(ret$fileID, "100kyna|100serine|ctrl|vehicle|APV"))) {
						ret$chemical_condition<- str_extract(ret$fileID, "100kyna|100serine|ctrl|vehicle|APV")	
						}
						
						###protocol call is very hacky
						#ret$timePoint<- str_extract(ret$fileID, "timept\\d")
						if (unique(str_detect(ret$fileID, "\\d\\dmin|\\dmin|t\\d\\d"))) {
							ret$timept_min <- str_extract(ret$fileID, "\\d\\dmin|\\dmin|t\\d\\d")
						}
						if (unique(str_detect(ret$fileID, "\\d\\d\\d\\ds"))) {
							ret$timept_s <- str_extract(ret$fileID, "\\d\\d\\d\\ds")
						}
						
						
						if (unique(str_detect(ret$fileID, "repl\\d\\d|repl\\d"))) {
							ret$exposeNum<-	str_extract(ret$fileID, "repl\\d\\d|repl\\d")
						}
						if (unique(str_detect(ret$fileID, "trx")) ) {
							ret$Trolox <- str_extract(ret$fileID, "notrx|trx")
							

						}	
						

						if (unique(str_detect(ret$fileID, "mMsuc")) ) {
							ret$sucrose<- str_extract(ret$fileID, "ctrl|\\d\\d\\dmMsuc|\\d\\dmMsuc")
							ret$sucrose_mM<- 	case_when( ret$sucrose == "ctrl" ~ 0,
																	ret$sucrose == "50mMsuc" ~ 50,
																	ret$sucrose == "100mMsuc" ~ 100,
																	ret$sucrose == "200mMsuc" ~ 200)
						}	
						ret$Ca<-	str_extract(ret$fileID, "\\d\\dmMCa|\\dmMCa|\\dpt\\d\\dmMCa|\\dpt\\dmMCa|\\dpt\\d\\dCa|\\dpt\\dCa|\\dCa")
													
																			

						ret$Ca_mM<-case_when(ret$Ca == "0pt25Ca" ~ 0.25,
																	ret$Ca == "0pt5Ca" ~ 0.5,
																	ret$Ca == "0pt6Ca" ~ 0.6,
																	ret$Ca == "2pt5Ca" ~2.5,
																	ret$Ca == "2Ca" ~ 2,
																	
																	ret$Ca == "8Ca" ~ 8,
																	ret$Ca == "0pt8Ca" ~ 0.8,
																	ret$Ca == "1pt2Ca" ~ 1.2,
																	ret$Ca == "1pt6Ca" ~ 1.6,
																	ret$Ca == "1Ca" ~ 1,
																	ret$Ca == "4Ca" ~ 4,

																	ret$Ca == "1pt2mMCa" ~ 1.2,
																	ret$Ca == "05mMCa" ~ 0.5,
																	ret$Ca == "1mMCa" ~ 1,
																	ret$Ca == "2mMCa" ~ 2,
																	ret$Ca == "4mMCa" ~ 4,
																	ret$Ca == "8mMCa" ~ 8
																	)


						#if (unique(str_detect(ret$fileID, "HNK")) ) {
						ret$HNK_phase<-str_extract(ret$fileID, "2pt5Ca_washout_long|2pt5Ca_washout|30uM-HNK_long_wait|30uM-HNK")
						#}

						ret$segmentation<- segmentationMethod

						if (unique(str_detect(ret$fileID, "blue\\d\\d_red\\d\\d")) ) {
							ret$bluepwr <- str_extract(ret$fileID, "blue\\d\\d")
							ret$redpwr <-  str_extract(ret$fileID, "red\\d\\d")
						}
						if (unique(str_detect(ret$fileID, "37C")) ) {
							ret$temperature <- str_extract(ret$fileID, "37C")
							
						} else {
							ret$temperature <-  "RT_but_check"
						}


						if (TTL != "nostim") {
							ret$TTL_start<- TTL

						}
						

						ret$ROINumber<-str_extract(ret$Source, "ROI\\d\\d\\d\\d")
						ret
		}


#ldply converts list object to data.frame
tic()
data<- ldply(filenames, read_csv_filename)

#select drops
drops <- c("Source","Ca_str",'timeStamp')  #extraneous data.frame columns clogging up memory
data<- data[ , !(names(data) %in% drops)]


if("Trolox" %in% colnames(data)) {
		data$Trolox[which(is.na(data$Trolox ) )] <- "no_trx"
	}
if("HNK_phase" %in% colnames(data)){
		data$HNK_phase[which(is.na(data$HNK_phase))] <- "ctrl"
}
data<- data %>% group_by_at(groupers) %>% mutate(interFrame = absoluteTime - lag(absoluteTime) )




print("Finished reading files into data.frame with necessary identifiers.")





toc()


