
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


#optional, if using script in isolation setwd within script 
#dir = setwd("\\\\blanpiedserver\\NASShare3\\Sam\\2022-02-07\\1Hz test")



### Changelog ###

### 5.23.22 ###

# commented out sucrose_str codes
# commented out DIV
# Re-integrated stimParadigm



## Read files into R 

filenames<- list.files(pattern=".csv", full.names=TRUE, recursive=TRUE) #Recurse through all subdirectories in the directory


#Define read_csv_filename, which names existing columns and adds new columns according to extracted strings
read_csv_filename <- function(filename){
		  ret <- read.csv(filename)
						colnames(ret)[1]<- "frameNumber"
					  colnames(ret)[2]<- "intensity"
						colnames(ret)[3]<- "fileID"
						colnames(ret)[4]<- "index"
						colnames(ret)[5]<- "timeStamp"
						colnames(ret)[6]<- "absoluteTime"
						
						
						ret$Source <- filename
						ret$sensor<- str_extract(ret$fileID, "jRGECO|GluSnFR3|GCaMP8f|GCaMP8m")
						ret$marker<-str_extract(ret$fileID, "nomarkers|SyPhy-mRuby3|PSD95FingR|Syn1a")								
						ret$dish<- str_extract(ret$fileID, "dish\\d")
						ret$plate<-str_extract(ret$fileID, "plate\\d\\d|plate\\d")
						ret$region<-str_extract(ret$fileID, "region\\d\\d|region\\d")
						

						ret$neuronSegment<- str_extract(ret$fileID, "axon|dendrite")
														
																				#			#				#			#		#
						ret$protocol<- 	str_extract(ret$fileID, "BasalRP|PP\\d\\d\\d\\d|PP\\d\\d\\d|PP\\d\\d|\\d\\dHz|\\dHz") 	
																	#str_detect( ret$fileID,"BasalRP | PP\\d\\d\\d\\d | PP\\d\\d\\d | PP\\d\\d|Train") ~  
																								#			#				#			#			#		#
																	#str_extract(ret$fileID, "BasalRP | PP\\d\\d\\d\\d | PP\\d\\d\\d | PP\\d\\d | \\d\\dHz | \\dHz"),
													#TRUE ~ "protocol not detected")

													#legacy
													#str_detect( ret$fileID, "PP\\d\\d\\d\\d | PP\\d\\d\\d | PP\\d\\d") ~ str_extract(ret$fileID, "PP\\d\\d\\d\\d | PP\\d\\d\\d | PP\\d\\d"),
													#str_detect( ret$fileID, "Train") ~ str_extract(ret$fileID, "\\d\\dHz | \\dHz"),
													
						
						###protocol call is very hacky

						ret$exposeNum<-	str_extract(ret$fileID, "_Glu_\\d|_Glu")
													

						ret$Ca<-	str_extract(ret$fileID, "\\d\\dmMCa|\\dmMCa")
													
																			

						ret$Ca_mM<-case_when( ret$Ca == "05mMCa" ~ 0.5,
																	ret$Ca == "1mMCa" ~ 1,
																	ret$Ca == "2mMCa" ~ 2,
																	ret$Ca == "4mMCa" ~ 4,
																	ret$Ca == "8mMCa" ~ 8)


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



muddle_detector<- function(df,df_colname_as_string) {
        
				findMuddle<- sum( as.integer( grepl(detector, df[[df_colname_as_string]]) ), na.rm=TRUE )
        
				
				if( findMuddle > 1) {
					print(paste0("Some datasets lacked strings for column: ", df_colname_as_string ) )
					
					
				} else {
				  print(paste0("The ...", df_colname_as_string, "... column did not have muddled observations." ) )
				}
				updatedMuddle = howMuddled + findMuddle
				assign("howMuddled", updatedMuddle, envir = .GlobalEnv)

}


print('Checking if any of the dataset is missing?')
detector <- "not detected"
howMuddled<- 0

muddle_detector(data, "neuronSegment")
muddle_detector(data, "protocol")
muddle_detector(data, "Ca")

if( howMuddled > 0) {
	totalObs = nrow(data)
	print(paste("The dataset is muddled. Muddle detector found", howMuddled, 
	            "observations, compared to the total number of observations,", totalObs, 
	            "... That's", round(howMuddled/totalObs*100, 2), "% of the current dataset!", sep = " "))

	#print("What would you like to do about the muddled entries? Type 0 to do nothing, or a string representing which column needs fixing.")
	#ui_getCol = readline()

	#	if(ui_getCol == "0") {
	#		print("Think on it!")
	#		break
	#		} else {
#
#				print("What's the plan? Type 0 to remove those observations or a string with which to replace muddled entries.")
#				ui_getEntry = readline()
#					if(ui_getEntry == "0") {
#							data_updated <- remove_obs(data, ui_getCol)
#						} else {
#							data_updated <- replace_obs(data,ui_getCol,ui_getEntry)
#					}
#				}
	} else {
	print("Looks like all the expected data is accounted for.")
}

print("Finished reading files into data.frame with necessary identifiers.")




rm(detector,howMuddled)

toc()




#########	Legacy code


#replace_obs<- function(df, var, replacement_string) {
#				data_updated<- df %>% dplyr::mutate(!! var := ifelse(grepl(detector, .), replacement_string, !!sym(var) ) )
#
#						!!sym(var) =  ifelse(str_detect(grepl(detector)) == TRUE, replacement_string, !!sym(var)  )    )
#				print("Replaced observations with: ", replacement_string)
#				data_updated

#}

#remove_obs<- function(df, var) {
#				data_updated<-df %>% dplyr::filter(!grepl(detector,!!sym(var)) ) 
#				print("Removed observations from dataframe.")
#				data_updated
#	}







#drops<- c("exposureTime", "laserPower", "paradigm", "stimFrequency", "Ca_mM") #3.20.22 producing figure for grant
#data<- data[ , !(names(data) %in% drops)]





#are.unique.colnames <- function(array){
#  return(length(unique(colnames(array))) == dim(array)[2])
#}

#ifelse(are.unique.colnames(data)==TRUE,
#				print("There are no duplicate colnames, this function may be useful in the future."),
#				print("WTF??!? There's a duplicate column?? That would literally break the dataframe.")
#				)




#should look to make this less hacky
			

			#ret$Source <- filename
			#ret$dirname <- dirname(filename) #may be unnecessary
			
			#Establish identifiers from column 3: fileID
			#ret$exposureTime<-str_extract(ret$fileID, "ET\\dms|ET\\d\\dms")
			#ret$laserPower<-str_extract(ret$fileID, "\\d\\d%|\\d%")
			
			#ret$testProtein<- str_extract(ret$Source, "PTP-cleavable|PTP-uncleavable|alone")
			#ret$Thrombin<- str_extract(ret$Source, "addThr|noThr|Thr")
			#ret$Thrombin_state<-case_when(ret$Thrombin == "addThr" ~ "+",
											#ret$Thrombin == "Thr" ~ "+",
											#ret$Thrombin == "noThr" ~ "-")

			#ret$marker<-str_extract(ret$fileID, "nomarkers|PSD95|RIMBP2|Syn1a")
			#ret$DIV<- str_extract(ret$fileID, "DIV\\d\\d")



			#ret$exposeNum<- case_when(ret$exposeFile == "_1_" ~ "1",
			#											ret$exposeFile == "_2_" ~ "2",
			#											ret$exposeFile == "_3_" ~ "3",
			#											ret$exposeFile == "_4_" ~ "4")
			#											ret$exposFile == "_5_" ~ "5",
			#											ret$exposFile == "_6_" ~ "6",
			#											ret$exposFile == "_7_" ~ "7",
			#											ret$exposFile == "_8_" ~ "8")
			
			#ret$stimFrequency<-str_extract(ret$fileID, "1Hz|5Hz")
			
		

										#str_extract( ret$fileID, "BasalRP |
										#			 PP\\d\\d\\d\\d | PP\\d\\d\\d | PP\\d\\d |
										#			 \\dHzTrain") ~ "")

										#ret$stim == "sham" ~ "Sham_stimulus",
										#	ret$stim == "pt5Hz" ~ "0.5_Hz",
										#	ret$stim == "1Hz" ~ "1_Hz",
										#	ret$stim == "2Hz" ~ "2_Hz",
										#	ret$stim == "5Hz" ~ "5_Hz")
			
			

			#for Ca panels
			
			#for sucrose panels
			#ret$sucrose<-str_extract(ret$fileID, "\\d\\d\\dmMsucrose|\\d\\dmMsucrose|\\dmMsucrose")
			#ret$sucrose_str<-case_when(ret$sucrose == "0mMsucrose" ~ "0 mM sucrose",
			#														ret$sucrose == "50mMsucrose" ~ "50 mM sucrose",
			#														ret$sucrose == "100mMsucrose" ~ "100 mM sucrose",
			#														ret$sucrose == "200mMsucrose" ~ "200 mM sucrose")
			#ret$sucrose_mM<-case_when(ret$sucrose == "0mMsucrose" ~ 0,
			#														ret$sucrose == "50mMsucrose" ~ 50,
			#														ret$sucrose == "100mMsucrose" ~ 100,
			#														ret$sucrose == "200mMsucrose" ~ 200)
