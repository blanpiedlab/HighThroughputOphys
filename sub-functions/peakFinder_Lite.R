#### Changelog 9.11.22 #### 
# 1. added a new param to plot in schmittTrig call - smoothed dFF. Use a span of 3 (high point density data from 9.8.22, 200 Hz sampling rate)



## peakFinder_Lite.R
## peakFinder_Lite.R will be essentially the same as peakFinder.R.
## peakFinder_Lite.R aims to take the generous peak indices (peak_idx) from peakFinder_Heavy.R to redefine the baseline with more precision. 
## By blending data.frame manipulations from read_and_flagOutliers.R, we will estimate a baseline using a window (25 or 50 frames?) where only
## the indices marked "notPeak" are used.  


#suppress warnings
oldw <- getOption("warn")
options(warn = -1)


library(gsignal)


dataset_id = unique(peaksHeavy$sensor)

if (length(dataset_id) > 1) {
		print("Are you sure you want to run this on multiple sensor datasets?")
		
		window = 100
		enter = 10
		exit = 20
		window_boundary = 25
		

	} else if (dataset_id == "GluSnFR3") {
		print(paste0("Flagging outliers for a dataset collected with a fast sensor : " dataset_id))
		
		window = 100
		enter = 10
		exit = 20
		window_boundary = 25
		
	} else if (dataset_id == "jRGECO|JF646|GCaMP8m") {
		print(paste0("Flagging outliers for a dataset collected with slow sensor : " dataset_id))
		
		window = 160
		enter = 20
		exit = 50
		window_boundary = 50
		

	} else if (dataset_id == "GCaMP8f") {
		print(paste0("Flagging outliers for a dataset collected with medium-speed sensor : " dataset_id))
		
		window = 100
		enter = 15
		exit = 35
		window_boundary = 35
		
	}




tic()

peaksLite <- peaksHeavy %>%
      	group_by_at(groupers) %>%
      	mutate(intensity_notPeak = ifelse(peak_idx == "notPeak"|flag == "notOutlier",intensity, NA),
      			baseline_v0 = na.locf(intensity_notPeak, fromLast = FALSE, na.rm = FALSE),
      			baseline = rollapply(baseline_v0, window, median, fill = NA, align="center")
      			)

								  
# Above code defines second-pass baseline adjustment. 
# Whenever an index is flagged "isPeak", intensity_notPeak returns NA. baseline_v0 back-fills these indices by replicating the last numeric value (step function at outliers)
# Finally, we smooth this using a rolling median. 


print("Finished estimating the baseline.")

toc()






tic()

## clean this up later
drops <- c("dFF","signal",'signal_median_v0')  #extraneous data.frame columns clogging up memory
peaksLite<- peaksLite[ , !(names(peaksLite) %in% drops)]



k.points = 3 #points to rollmean 

dFF <- peaksLite %>%
      	group_by_at(groupers) %>%
      	mutate(normIntensity = intensity/max(intensity,na.rm=TRUE),
      				normBaseline = baseline/max(intensity,na.rm=TRUE),
      				F = intensity,
      				dF = (intensity - baseline),
      				dFF = (intensity - baseline)/baseline,
      				smooth_dFF = rollapply(dFF, k.points, mean, fill = NA, align="center"),
      				
      				dFF_idx = as.numeric(ifelse(peak_idx=="notPeak", dFF, "NA")),
    				dFF_signal_v0 = na.locf(dFF_idx, fromLast = TRUE, na.rm = FALSE),				
      				outlier_indices = case_when(peak_idx == "isPeak" ~ "blue", 
      											flag == "Outlier" ~ "red") 
      				) 


getTrace.stdev = suppressMessages(dFF %>% group_by_at(groupers) %>%
					dplyr::filter(peak_idx == "notPeak" & flag == "notOutlier") %>%
					summarise(dFF_stdev = sd(dFF_idx,na.rm=TRUE))
					)

getTrace.yLimits = suppressMessages(dFF %>% summarise(max_dFF = max(dFF,na.rm=TRUE),
										min_dFF = min(dFF, na.rm=TRUE))
									)


dFF <- suppressMessages(left_join(dFF,getTrace.stdev))
dFF<- suppressMessages(left_join(dFF,getTrace.yLimits))
drops <- c('dFF_idx',"dFF_signal_v0")  #extraneous data.frame columns clogging up memory
dFF<- dFF[ , !(names(dFF) %in% drops)]


print("Finished calculating dFF.")

toc()



tic()
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v4 - markerSegmentation/schmittTrig_v4.R") ## import Schmitt Trigger function with plotting









#from Stack Overflow: https://stackoverflow.com/questions/5135087/change-value-of-some-column-in-xts-based-on-other-columns-values-with-lookback
# Defines any positive cumulative sum as a 1. Useful for defining a region using some start and stop signal as the boundaries.

cumplus <- function(y) Reduce(function(a,b) a + b > 0, y, 0, accum=TRUE)[-1]


# threshold over which to identify a peak. Lower trigger = 1.5*sigma, upper trigger = 3.5*sigma
threshold = c(1.5,5)    
thresh_stdv = 1.5

peaks<-dFF %>%
		group_by_at(groupers) %>%
			do(schmittTrig_v4(dataframe=., 
								time=.$absoluteTime, 
								dFF=.$dFF, 
								normIntensity=.$normIntensity, 
								normBaseline=.$normBaseline, 
								outlier_indices=.$outlier_indices,
								threshold=threshold, 
								std=.$dFF_stdev,
								ymax=.$max_dFF,
								ymin=.$min_dFF, 
								thresh_stdev=thresh_stdv,
								smoothed_dFF = .$smooth_dFF)
					) %>% 	  
			mutate(interFrame = absoluteTime - lag(absoluteTime),
					ShortEntry = case_when(signal == 1 ~ 0,
											lead(signal, n= enter) == 1 ~ 1,
											signal == 0 ~ 0),
					ShortExit = case_when(signal == 1 ~ 0,
											lag(signal, n= exit) == 1 ~ 1,
											signal == 0 ~ 0),
					cumplus = cumplus(ShortEntry - ShortExit),
					temp = cumplus - c(0,pmin(0,diff(cumplus))),
					peakIncrement = cumsum(c(0, as.numeric(diff(temp)>=1))),
					temp2 = case_when(temp == 1 ~ peakIncrement,
										temp == 0 ~ 0),
					peakID = case_when(temp2 == 0 ~ "NotPeak",
										temp2 > 0 ~ paste0("Peak",temp2)),
					WindowEntry = case_when(ShortEntry >= 1 ~ 0,
											lead(ShortEntry, n= window_boundary) >= 1 ~ 1,
											ShortEntry == 0 ~ 0),
					WindowExit = case_when(ShortExit >= 1 ~ 0,
											lag(ShortExit, n= window_boundary) >= 1 ~ 1,
											ShortExit == 0 ~ 0),
					cumplusWindow = cumplus(WindowEntry - WindowExit),
					temp3 = cumplusWindow - c(0,pmin(0,diff(cumplusWindow))),
					windowIncrement = cumsum(c(0, as.numeric(diff(temp3)>=1))),
					windowedPeakID = case_when(temp3 == 0 ~ "NotPeak",
										temp3 > 0 ~ paste0("windowedPeak",windowIncrement))

						) %>% 
			select(-temp, -temp2, -peakIncrement) #get rid of extraneous columns



					# technically don't need to peakIncrement in peaksHeavy. Should troubleshoot later.
					
			
drops <- c("temp","temp2","peakIncrement", "normIntensity", "peak_idx","intensity_notPeak","baseline","outlier_indices","ShortExit","ShortEntry","cumplus","normStDev","normBaseline", "baseline_v0", "F", "dF", "dFF_stdev","signal_stdev","flag","frameNumber",
			"WindowEntry","WindowExit","cumplusWindow","temp3")  #extraneous data.frame columns clogging up memory
peaks<- peaks[ , !(names(peaks) %in% drops)] %>%
		na.omit()

#warnings back on
options(warn = oldw)

rm(list=ls()[! ls() %in% c("peaks","groupers")])  #clear all vars in memory except for flagged data.

										 

print("Finished identifying peaks and assigning peakIDs for analysis.")
toc()



print("Finished running the data extraction portion of pipeline v5.0.")
