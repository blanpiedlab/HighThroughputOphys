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
interFrame = mean(peaksHeavy$interFrame,na.rm=TRUE)
minimum_window_duration = 1.0
k=3


if (length(dataset_id) > 1) {
		print("Are you sure you want to run this on multiple sensor datasets?")
		
		if(any(dataset_id %in% c("JF646", "jRGECO")) ) {
					rise_time = 0.2
					full_decay_time = 2.0
					wb_val = 0.50

					window = ceiling(minimum_window_duration/interFrame)
					enter = ceiling(2*rise_time/interFrame) 
					exit = ceiling(2*full_decay_time/interFrame) 
					window_boundary = ceiling(wb_val/interFrame)
				} else {
					rise_time = 0.0193
					full_decay_time = 0.08
					wb_val = 0.1

					window = ceiling(minimum_window_duration/interFrame)
					enter = ceiling(2*rise_time/interFrame) 
					exit = ceiling(2*full_decay_time/interFrame) 
					window_boundary = ceiling(wb_val/interFrame)



				}

	} else if (dataset_id == "GluSnFR3") {
		print(paste0("Flagging outliers for a dataset collected with a fast sensor : ", dataset_id))
		rise_time = 0.0193
		full_decay_time = 0.08
		wb_val = 0.1


		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2*rise_time/interFrame) 
		exit = ceiling(2*full_decay_time/interFrame) 
		window_boundary = ceiling(wb_val/interFrame)

	} else if (dataset_id == "jRGECO" | dataset_id == "JF646" | dataset_id == "GCaMP8m") {
		print(paste0("Flagging outliers for a dataset collected with slow sensor : ", dataset_id))
		rise_time = 0.075
		full_decay_time = 0.25
		wb_val = 0.50

		
		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2*rise_time/interFrame) 
		exit = ceiling(2*full_decay_time/interFrame) 
		window_boundary = ceiling(wb_val/interFrame)


	} else if (dataset_id == "GCaMP8f") {
		print(paste0("Flagging outliers for a dataset collected with medium-speed sensor : ", dataset_id))
		rise_time = 0.025
		full_decay_time = 0.150
		wb_val = 0.1

		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2*rise_time/interFrame) 
		exit = ceiling(2*full_decay_time/interFrame) 
		window_boundary = ceiling(wb_val/interFrame)


	}




tic()

peaksLite <- peaksHeavy %>%
      	group_by_at(groupers) %>%
      	mutate(intensity_notPeak = as.numeric(ifelse(peak_idx == "notPeak",intensity_smooth, NA)),
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
      	mutate(normIntensity = intensity_smooth/max(intensity_smooth,na.rm=TRUE),
      				normBaseline = baseline/max(intensity_smooth,na.rm=TRUE),
      				F = intensity_smooth,
      				dF = (intensity_smooth - baseline),
      				dFF = (intensity_smooth - baseline)/baseline,
      				smooth_dFF = rollapply(dFF, k.points, mean, fill = NA, align="center"),
      				
      				dFF_idx = as.numeric(ifelse(peak_idx=="notPeak", dFF, "NA")),
    					dFF_signal_v0 = na.locf(dFF_idx, fromLast = FALSE, na.rm = FALSE),				
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
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/schmittTrig_v4.R") ## import Schmitt Trigger function with plotting









#from Stack Overflow: https://stackoverflow.com/questions/5135087/change-value-of-some-column-in-xts-based-on-other-columns-values-with-lookback
# Defines any positive cumulative sum as a 1. Useful for defining a region using some start and stop signal as the boundaries.

cumplus <- function(y) Reduce(function(a,b) a + b > 0, y, 0, accum=TRUE)[-1]


# threshold over which to identify a peak. Lower trigger = 2*sigma, upper trigger = 3.5*sigma
threshold = c(1.5,5)    
thresh_stdv = 1.5

peaks<-dFF %>%
		group_by_at(groupers) %>%
			do(schmittTrig_v4(dataframe=., 
								time=.$index, 
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
					




if(savePlots == TRUE) {

	library(ggeasy)
	library(scales)

	subDir <- paste0("peakFinderOutputs",first(unique(peaks$sensor) ) )
  


	#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme_peakFinder_outputs.R") 
	

	#do things
	#want to save out schmitt Trigg 3 essentially

	.ROI_keyvars = rlang::syms(ROI_keys)
	tmp_df<- peaks %>% mutate(ROI_keys = paste(!!!.ROI_keyvars, sep="-"))

	ROIs = unique(tmp_df$ROI_keys)
	ROIs_to_sample = ceiling( length(ROIs) /n ) 
	randomROIs<- sample(ROIs,ROIs_to_sample)

	print( paste0("There are ", length(ROIs), " total unique ROIs in this dataset. We are going to save out traces for ", (1/n)*100, "% of them.") )

	 for (i in 1:length(randomROIs)) {

		tmp_data <- tmp_df %>%
        dplyr::filter(ROI_keys == randomROIs[i]) 
        lvl = mean(tmp_data$dFF_stdev, na.rm = TRUE)*threshold
    	ymax = mean(tmp_data$max_dFF, na.rm=TRUE)
    	ymin = mean(tmp_data$min_dFF, na.rm=TRUE)
   	






   
	## p1 - raw trace wwith outliers
	## p2 - normalized intensity with baseline overlay
	## p3 - dFF calculated with schmitTrig levels as geom_hline
	## p4 - signal binary

	p1<-ggplot(tmp_data, aes(x=index, y = normIntensity)) +
            geom_line(aes(fill='Normalized Signal'),size = 3) +
            geom_point(aes(fill="Outliers Detected",colour=outlier_indices), size=5,shape=21,fill="white",stroke=2,alpha=0.5)+
            labs( x="Index",
                  y="Normalized Intensity",
                  title=unique(tmp_data$ROI_keys),
                  subtitle="Normalized Trace with Outliers Identified")+
            theme_tufte()+
            scale_colour_identity(guide="legend")+
            scale_y_continuous(labels = number_format(accuracy = 0.1))+
            #coord_cartesian(ylim=c(ymin,ymax))+
            #scale_colour_manual()+
            scale_fill_manual(values=c('Normalized Signal'='black','Outliers Detected' = NA))+
            my.theme+
            theme(plot.title = element_text(size = 15),
            		legend.position="none")+
            easy_remove_x_axis()


    p2<-ggplot(tmp_data, aes(x=index, y = normIntensity)) +
            geom_line(aes(colour="Normalized Signal"),size = 3 ) +
            geom_line(aes(x=index,y=normBaseline,colour="Estimated Baseline"),size=3.5) + 
            
            labs( x="Index",
                  y="Normalized Intensity",
                  title="",
                  subtitle="Normalized Trace with Estimated Baseline")+
            theme_tufte()+
            scale_y_continuous(labels = number_format(accuracy = 0.1))+
            
            scale_colour_manual(values=c('Normalized Signal'='black','Estimated Baseline'='red'))+
            my.theme+
            theme(legend.position="none")+
            easy_remove_x_axis()   


    p3<-ggplot(tmp_data, aes(x=index, y = dFF)) +
            geom_line(aes(colour="Normalized Signal"),size = 3 ) +
            geom_hline(colour="blue",yintercept=lvl,size=2) + 
            
            labs( x="Index",
                  y=expression(Delta*"F/F"),
                  title="",
                  subtitle=expression(Delta*"F/F with Schmitt Trigger"))+
            theme_tufte()+
            coord_cartesian(ylim=c(ymin,ymax))+
            scale_y_continuous(labels = number_format(accuracy = 0.1))+
            scale_colour_manual(values=c('Normalized Signal'='black',"Schmitt Trigger" ='blue'))+
            my.theme+
            theme(legend.position="none")+
            easy_remove_x_axis()    

     p4<-ggplot(tmp_data, aes(x=index, y = signal)) +
             geom_line(aes(colour="Extracted Signal Indices"),size = 3 ) +
            
             labs( x="Index",
                   y="Signal",
                   title="",
                   subtitle="")+
             theme_tufte()+
             scale_y_continuous(breaks=c(0.0,1.0),labels = number_format(accuracy = 0.1))+
             coord_cartesian(ylim=c(-0.1,1.1))+

            
             scale_colour_manual(values=c('Extracted Signal Indices'='red'))+
             my.theme+
             theme(legend.position="none",
             		plot.title=element_blank())

     fullplot=grid.arrange(p1, p2, p3,p4, nrow=4)



	ggsave(filename=paste0("peakFinder_trace#",i,".png"),plot=fullplot, device="png",dpi=600, units="in",width=12,height=20)
   



rm(tmp_data)
}
      

}
rm(tmp_df)
			
drops <- c("temp","temp2","peakIncrement", "normIntensity", "peak_idx","intensity_notPeak","baseline","outlier_indices","ShortExit","ShortEntry","cumplus","normStDev","normBaseline", "baseline_v0", "F", "dF", "dFF_stdev","signal_stdev","flag","frameNumber",
			"WindowEntry","WindowExit","cumplusWindow","temp3")  #extraneous data.frame columns clogging up memory
peaks<- peaks[ , !(names(peaks) %in% drops)] %>%
		na.omit()

#warnings back on
options(warn = oldw)

rm(list=ls()[! ls() %in% c("peaks")])  #clear all vars in memory except for flagged data.

										 

print("Finished identifying peaks and assigning peakIDs for analysis.")
toc()



print("Finished running the data extraction portion of pipeline v5.0.")
