
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
						colnames(ret)[7]<- "sensor"
						
						#ret$Source <- filename
						
						
						#GluSnFR3_dish1_plate01_region5_repl01_ROI0001
						ret$fileID <- filename
						#ret$infectionID_Glu <-str_extract(ret$fileID, "\\d\\dkGlu|\\dkGlu")
						#ret$infectionID_spine <-str_extract(ret$fileID, "\\d\\dkJF|\\dkJF")
						#ret$chemical_condition<- str_extract(ret$fileID, "ctrl|APV|\\d\\d\\duMMg|\\d\\duMMg")
						#ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|JF646|GCaMP8f|GCaMP8m")
						ret$dish<-str_extract(ret$fileID, "dish\\d\\d|dish\\d")
						ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")							
						ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
						ret$replicate<-str_extract(ret$fileID, "repl\\d\\d|repl\\d")
						ret$ROINumber<- str_extract(ret$fileID, "ROI\\d\\d\\d\\d")
						ret
		}


#ldply converts list object to data.frame
tic()
data<- ldply(filenames, read_csv_filename)

#select drops
drops <- c("Source","Ca_str",'timeStamp')  #extraneous data.frame columns clogging up memory
data<- data[ , !(names(data) %in% drops)]
data<- data %>% group_by_at(groupers) %>% mutate(interFrame = absoluteTime - lag(absoluteTime) )


print("Finished reading files into data.frame with necessary identifiers.")





toc()
