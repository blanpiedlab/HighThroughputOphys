
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
						
						
						ret$Source <- filename
						ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|JF646|GCaMP8f|GCaMP8m")
						ret$transduction <- str_extract(ret$fileID, "TRANSFECT|LENTI")
						#ret$marker<-str_extract(ret$fileID, "nomarkers|SyPhy-mRuby3|PSD95FingR|Syn1a")								
						ret$dish<- str_extract(ret$fileID, "dish\\d")
						ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")
						ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
						

						ret$neuronSegment<- str_extract(ret$fileID, "axon|dendrite")
														
																				#			#				#			#		#
						ret$protocol<- 	str_extract(ret$fileID, "BasalRP|PP\\d\\d\\d\\d|PP\\d\\d\\d|PP\\d\\d|\\d\\dHz|\\dHz") 	
						
						ret$exposeNum<-	str_extract(ret$fileID, "repl\\d|_Glu")
													

						ret$Ca<-	str_extract(ret$fileID, "\\d\\dmMCa|\\dmMCa")
													
																			

						ret$Ca_mM<-case_when( ret$Ca == "05mMCa" ~ 0.5,
																	ret$Ca == "1mMCa" ~ 1,
																	ret$Ca == "2mMCa" ~ 2,
																	ret$Ca == "4mMCa" ~ 4,
																	ret$Ca == "8mMCa" ~ 8)
						ret$segmentation<- segmentationMethod
						ret$TTL_start<- TTL


						ret$ROINumber<-str_extract(ret$Source, "ROI\\d\\d\\d\\d")
						ret
		}


#ldply converts list object to data.frame
tic()
data<- ldply(filenames, read_csv_filename)

#select drops
drops <- c("Source","Ca_str",'timeStamp')  #extraneous data.frame columns clogging up memory
data<- data[ , !(names(data) %in% drops)]

data$exposeNum[which(data$exposeNum == "_Glu")]<- "_Glu_0"
#data$marker[which(data$marker == NA)]<- "nomarkers"
data<- data %>% group_by_at(groupers) %>% mutate(interFrame = absoluteTime - lag(absoluteTime) )



print("Finished reading files into data.frame with necessary identifiers.")






toc()



