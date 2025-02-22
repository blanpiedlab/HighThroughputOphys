##flagOutliers.R


#suppress warnings
oldw <- getOption("warn")
options(warn = -1)

# This code defines first-pass baseline adjustment. 
# Whenever an outlier is flagged, baseline_idx returns NA. baseline_v0 back-fills these indices by replicating the last numeric value (step function at outliers)
# This results in the pseudo-baseline from which we will estimate dFF. Following peak identification, we will refine the baseline index estimation. 



# this section of code defines window in a generalized fashion



dataset_id = unique(data$sensor)
interFrame = mean(data$interFrame,na.rm=TRUE)
minimum_window_duration = 0.75 #seconds
k_duration_normal = 0.04 #seconds
k_duration_long = 0.04 #seconds


if (length(dataset_id) > 1) {
	print("Are you sure you want to run this on multiple sensor datasets?")
	k=ceiling(k_duration_normal/interFrame)
	window = ceiling(minimum_window_duration/interFrame)
} else if (dataset_id == "JF646"){
	print(paste0("Flagging outliers for a dataset collected with sensor : ", dataset_id) )
	k=1#ceiling(k_duration_long/interFrame)
	window = ceiling(minimum_window_duration/interFrame)
} else if (dataset_id == "GluSnFR3") {
	print(paste0("Flagging outliers for a dataset collected with sensor : ", dataset_id) )
	k=1#ceiling(k_duration_normal/interFrame)
	window = ceiling(minimum_window_duration/interFrame)

} else {
	print(paste0("Flagging outliers for a dataset collected with sensor : ", dataset_id) )
	k=ceiling(k_duration_normal/interFrame)
	window = ceiling(minimum_window_duration/interFrame)

}
thresh_stdv = 1.0 #threshold over which to flag as outlier (this is a very generous threshold - captures most of the signal amplitude such that these aren't figured into the pseudo-baseline calculation)




data_smooth <- suppressMessages(data %>%
      group_by_at(groupers) %>%
      mutate(intensity_smooth = intensity,#rollapply(intensity, k, mean, fill=NA, align="center"),
      		signal_median_v0=rollapply(intensity_smooth, window, median, fill=NA, align="center")) )


getTrace.stdev = suppressMessages(data_smooth %>% group_by_at(groupers) %>%
					#group_by(laserPower, exposureTime, plate, region, Ca_mM, exposeNum, ROINumber) %>%
					summarise(signal_stdev = sd(intensity_smooth, na.rm=TRUE))
					)

data_flagged = suppressMessages(left_join(data_smooth,getTrace.stdev) %>% 

			group_by_at(groupers) %>%
			mutate(
      		#,
									
			#outlier classification - outlier anytime signalSmooth - signal_median > 1.5*stdev of the trace
			flag=ifelse((intensity_smooth - signal_median_v0) < (thresh_stdv * signal_stdev), "notOutlier","Outlier"), 
			
			#define pseudo_baseline
			pseudoBaseline_idx = as.numeric(ifelse(flag=="notOutlier", intensity_smooth, NA)),
			pseudoBaseline_v0 = na.locf(pseudoBaseline_idx, fromLast = FALSE, na.rm = FALSE),
			pseudoBaseline = rollapply(pseudoBaseline_v0, window, median, fill=NA, align="center"),
			)
		)  




### UNCOMMENT THIS BLOCK OF CODE TO SAVE PLOTS FROM flagOutliers.R ######




# if(savePlots == TRUE) {

# 	library(ggeasy)
# 	library(scales)

# 	if(is.null(data_flagged$protocol)) {
# 		subDir <- paste0("peakFinder_flagOutliers_outputs")	
# 	} else {
# 		subDir <- paste0("peakFinder_flagOutliers_outputs",first(unique(data_flagged$protocol) ) )
  	
# 	}
	


# 	#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
# 	source("YourDrive:\\HighThroughputOphys-main/sub-functions/setFigureDirectory.R") 
# 	source("YourDrive:\\HighThroughputOphys-main/sub-functions/myTheme_peakFinder_outputs.R") 
	

# 	#do things
# 	#want to save out schmitt Trigg 3 essentially

# 	.ROI_keyvars = rlang::syms(ROI_keys)
# 	tmp_df<- data_flagged %>% mutate(ROI_keys = paste(!!!.ROI_keyvars, sep="-"))

		
	
# 	ROIs = unique(tmp_df$ROI_keys)
# 	ROIs_to_sample = ceiling( length(ROIs) /n ) 
# 	randomROIs<- sample(ROIs,ROIs_to_sample)

	
# 	print( paste0("There are ", length(ROIs), " total unique ROIs in this dataset. We are going to save out traces for ", (1/n)*100, "% of them.") )

# 	 for (i in 1:length(randomROIs)) {

# 		tmp_data <- tmp_df %>%
#         			dplyr::filter(ROI_keys == randomROIs[i]) 
#         # lvl = mean(tmp_data$dFF_stdev, na.rm = TRUE)*threshold
# 	  #   	ymax = mean(tmp_data$max_dFF, na.rm=TRUE)
#     	# 	ymin = mean(tmp_data$min_dFF, na.rm=TRUE)
   	






   
# 	## p1 - raw trace with stdev threshold
# 	## p2 - normalized intensity with baseline overlay
# 	## p3 - dFF calculated with schmitTrig levels as geom_hline
# 	## p4 - signal binary

# 	p1<-ggplot(tmp_data, aes(x=absoluteTime, y = intensity_smooth, colour="black")) +
# 		geom_line(size=3.5)+
# 		geom_line(aes(x=absoluteTime,y=signal_median_v0,colour="red"),size=3.5) + 
#             geom_line(aes(x=absoluteTime,y=signal_median_v0+(thresh_stdv * signal_stdev),colour="green"),size=3.5) + 
#             #geom_line(aes(x=absoluteTime,y=signal_median_v0-(thresh_stdv * signal_stdev),colour="green"),size=3.5) + 
            
#             labs( x="Time (s)",
#                   y="Intensity",
#                   title=unique(tmp_data$ROI_keys),
#                   subtitle="Trace with threshold identification")+
#             theme_tufte()+
#             scale_colour_identity(guide="legend")+
#             scale_y_continuous(labels = number_format(accuracy = 0.1))+
#             #coord_cartesian(ylim=c(ymin,ymax))+
#             #scale_colour_manual()+
#             #scale_fill_manual(values=c('Normalized Signal'='black','Outliers Detected' = NA))+
#             my.theme+
#             theme(plot.title = element_text(size = 15),
#             		legend.position="none")#+
            

# 	ggsave(filename=paste0("peakFinder_flagOutliers_trace#",i,".png"),plot=p1, device="png",dpi=600, units="in",width=12,height=6)
   

# # p2 recolor things based on whether it should be considered outlier?



# # fullplot=grid.arrange(p1, p2, nrow=2)


# rm(tmp_data)
# }
      

# }
# rm(tmp_df)








drops <- c("pseudoBaseline_idx","pseudoBaseline_v0")  #extraneous data.frame columns clogging up memory
data_flagged<- data_flagged[ , !(names(data_flagged) %in% drops)]

								  
print("Finished flagging outliers and establishing a pseudo-baseline.")

toc()



#warnings back on
options(warn = oldw)


rm(list=ls()[! ls() %in% c("data_flagged","groupers",'showGraphs','mainDir','path','n','savePlots','randomROIs','ROIs','ROIs_to_sample','toggle_short_peaks')])  #clear all vars in memory except for flagged data.
if(!is.null(dev.list())) dev.off()

