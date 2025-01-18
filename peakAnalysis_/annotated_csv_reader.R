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



#Define read_csv_filename, which names existing columns and adds new columns according to extracted strings
annotated_csv_reader <- function(filename){
		  ret <- read.csv(filename)
						colnames(ret)[1]<- "rowNumber"
					  	ret
		}


#ldply converts list object to data.frame
tic()
data<- ldply(filenames, annotated_csv_reader)

