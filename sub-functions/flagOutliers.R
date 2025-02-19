##flagOutliers.R


#suppress warnings
oldw <- getOption("warn")
options(warn = -1)

# This code defines first-pass baseline adjustment. 
# Whenever an outlier is flagged, baseline_idx returns NA. baseline_v0 back-fills these indices by replicating the last numeric value (step function at outliers)
# This results in the pseudo-baseline from which we will estimate dFF. Following peak identification, we will refine the baseline index estimation. 



# this section of code defines window in a generalized fashion



dataset_id = unique(data$sensor)

if (length(dataset_id) > 1) {
	print("Are you sure you want to run this on multiple sensor datasets?")
	window = 200
} else {
	print(paste0("Flagging outliers for a dataset collected with sensor : " dataset_id))
	window = 200
}

thresh_stdv = 1.5 #threshold over which to flag as outlier (this is a very generous threshold - captures most of the signal amplitude such that these aren't figured into the pseudo-baseline calculation)



#collapse stdev to single value for every trace. Prevents stdev from dispersing wildly when large peaks occur. Makes outlier estimation cleaner and hopefully reduces dFF wobble  
getTrace.stdev = suppressMessages(data %>% group_by_at(groupers) %>%
					#group_by(laserPower, exposureTime, plate, region, Ca_mM, exposeNum, ROINumber) %>%
					summarise(signal_stdev = sd(intensity))
					)

tic()

data_flagged <- suppressMessages(left_join(data,getTrace.stdev) %>%
      group_by_at(groupers) %>%
      mutate(signal_median_v0=rollapply(intensity, window, median, fill=NA, align="center"),
									
			#outlier classification - outlier anytime signalSmooth - signal_median > 1.5*stdev of the trace
			flag=ifelse((abs(intensity - signal_median_v0) < (thresh_stdv * signal_stdev)), "notOutlier","Outlier"), 
			
			#define pseudo_baseline
			pseudoBaseline_idx = as.numeric(ifelse(flag=="notOutlier", intensity, NA)),
			pseudoBaseline_v0 = na.locf(pseudoBaseline_idx, fromLast = FALSE, na.rm = FALSE),
			pseudoBaseline = rollapply(pseudoBaseline_v0, window, median, fill=NA, align="center"),
			)
		)  

drops <- c("pseudoBaseline_idx","pseudoBaseline_v0")  #extraneous data.frame columns clogging up memory
data_flagged<- data_flagged[ , !(names(data_flagged) %in% drops)]

								  
print("Finished flagging outliers and establishing a pseudo-baseline.")

toc()

#warnings back on
options(warn = oldw)


rm(list=ls()[! ls() %in% c("data_flagged","groupers",'showGraphs')])  #clear all vars in memory except for flagged data.


