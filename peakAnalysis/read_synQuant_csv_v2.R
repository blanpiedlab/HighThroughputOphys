

# read_synQyantResults_csv.R

# read_synQuantResults_csv.R is a script to read in the SynQuant files that are produced from the marker-based segmentation imageJ macro. 
# the Ch2Results.csv files are in a deep directory within : Y:\Sam\2022-09-08\imageFiles\markerBased_analysis\zstackOutputs\


# the variables are in the files as follows: 
# ROINumber as an index,	Area,	Mean,	X,	Y,	Circ.,	IntDen,	RawIntDen,	AR,	Round,	Solidity,	Image Number,	Image Name

# something got borked 10/25/22


# if spont
# 	 	Min	Max	X	Y	IntDen	RawIntDen	Image Number	Image Name



#if singleAP
# Area	Min	Max	X	Y	IntDen	RawIntDen	Slice	Image Number	Image Name


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



SQ_reader <- function(synQuantDir, file_pattern,protocol){

		Ca_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )  
						

		setwd(synQuantDir)


		SQfilenames<- list.files(pattern=file_pattern, full.names=TRUE, recursive=TRUE) #Recurse through all subdirectories in the directory

		#print(paste0("This is the list of filenames we found :", SQfilenames))

		colname_switch = str_extract(synQuantDir, "singleAP|spont")

			#

			#Define read_csv_filename, which names existing columns and adds new columns according to extracted strings
			read_csv_filename <- function(filename, col_switch=colname_switch){
					  					  	
					  if(col_switch == "singleAP"){
					  						ret <- read.csv(filename)
					  						colnames(ret)[1]<- "ROINumber_asIndex"
											colnames(ret)[2]<- "Area"
											colnames(ret)[3]<- "Min"
											colnames(ret)[4]<- "Max"
											colnames(ret)[5]<- "X"
											colnames(ret)[6]<- "Y"
											colnames(ret)[7]<- "IntDen"	
											colnames(ret)[8]<- "RawIntDen"
											colnames(ret)[9]<- "Slice"
											colnames(ret)[10]<- "imageNumber"
											colnames(ret)[11]<- "fileID"
											
											
											ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|JF646|GCaMP8f|GCaMP8m")
											#ret$marker<-str_extract(ret$fileID, "SyPhy")
											ret$segmentation<-"marker"								
											ret$dish<- str_extract(ret$fileID, "dish\\d\\d|dish\\d")
											ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")
											ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
											#ret$Ca<-	str_extract(ret$fileID, "\\dpt\\dCa|\\dCa")
											ret$protocol_fixed<- protocol
																		
											ret$ROINumber<-sprintf("ROI%04d",ret$ROINumber_asIndex)
											ret

						 } else if(col_switch == "spont"){
						 					ret <- read.csv(filename)
					  						colnames(ret)[1]<- "ROINumber_asIndex"
						 					colnames(ret)[2]<- "Min"
											colnames(ret)[3]<- "Max"
											colnames(ret)[4]<- "X"
											colnames(ret)[5]<- "Y"
											colnames(ret)[6]<- "IntDen"
											colnames(ret)[7]<- "RawIntDen"
											colnames(ret)[8]<- "imageNumber"
											colnames(ret)[9]<- "fileID"
											
											
											ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|JF646|GCaMP8f|GCaMP8m")
											#ret$marker<-str_extract(ret$fileID, "SyPhy")
											ret$segmentation<-"marker"								
											ret$dish<- str_extract(ret$fileID, "dish\\d\\d|dish\\d")
											ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")
											ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
											ret$Ca<-	as.character(str_extract(ret$fileID, "\\dpt\\dCa|\\dCa") ) 
											ret$Ca<- factor(ret$Ca,levels = c("0pt5Ca", "1Ca", "2Ca", "4Ca"))
											levels(ret$Ca) <- Ca_levels
											
											ret$protocol_fixed<- protocol
																		
											ret$ROINumber<-sprintf("ROI%04d",ret$ROINumber_asIndex)


											ret



						 }
					}


			#ldply converts list object to data.frame
			tic()
			SQ_data<- ldply(SQfilenames, read_csv_filename)


			drop_cols<- c("ROINumber_asIndex","imageNumber","fileID","Slice","Area")
			SQ_data<- SQ_data[,!(names(SQ_data) %in% drop_cols)]                    
			


			print("Finished reading in the SynQuant data.")

			
			toc()
			SQ_data

}


