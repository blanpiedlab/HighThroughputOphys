## peakFinder-Heavy.R
 


# 	At the end of writing this script, I should consider at which breakpoints which columns can be dropped. (e.g. flag, pseudoBaseline, etc)

## The dFF for every ROI will be calculated using the pseudo-baseline in data_flagged (pseudoBaseline).  
## Then, peak identification is performed. To enable the user to watch the progress of the script, 
## schmittTrig.R is called, which flags signals and writes 4 plots to the cache. 
# Plot 1. Normalized signal intensity with normalized pseudoBaseline overlay. 
# Plot 2. Close-up of normalized pseudoBaseline.
# Plot 3. dFF with Schmitt Trigger threshold overlay. 
# Plot 4. Binary state in time (schmitt trigger signal ID)  

## The output from this script should be a data.frame with flexibly defined signal epochs. 
## The traceState column will provide a reference for further baseline adjustment in peakFinder-Lite.   



#suppress warnings
oldw <- getOption("warn")
options(warn = -1)

library(gsignal) ### time series analysis, necessary for Schmitt trigger identification of peaks

tic()


dataset_id = unique(data_flagged$sensor)
if(unique(data_flagged$segmentation) %in% c("activityGlu","marker")){
	print('toggling short peaks')
	toggle_short_peaks = TRUE
} else {toggle_short_peaks = FALSE }

interFrame = mean(data_flagged$interFrame,na.rm=TRUE)
minimum_window_duration = 0.5
k=3

if (length(dataset_id) > 1) {
		print("Are you sure you want to run this on multiple sensor datasets?")
		
		if(any(dataset_id %in% c("JF646", "jRGECO","GCaMP8m")) ) {
						rise_time = 0.2
						full_decay_time = 1.25
		
						window = ceiling(minimum_window_duration/interFrame)
						enter = ceiling(2.5*rise_time/interFrame)
						exit = ceiling(2.5*full_decay_time/interFrame)
				} else {
						rise_time = 0.0193
						full_decay_time = 0.10

						window = ceiling(minimum_window_duration/interFrame)
						enter = ceiling(2.5*rise_time/interFrame)
						exit = ceiling(2.5*full_decay_time/interFrame)



				}


	} else if (dataset_id == "GluSnFR3") {
		print(paste0("Flagging outliers for a dataset collected with a fast sensor : ", dataset_id))
		rise_time = 0.0193
		full_decay_time = 0.120


		if(toggle_short_peaks == TRUE){
					window = ceiling(minimum_window_duration/interFrame)
					enter = ceiling(0.5*rise_time/interFrame) 
					exit = ceiling(0.5*full_decay_time/interFrame) 
					
		} else {
		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2.5*rise_time/interFrame) 
		exit = ceiling(2.5*full_decay_time/interFrame) 
		}


		
	} else if (dataset_id == "jRGECO" | dataset_id == "JF646" | dataset_id == "GCaMP8m") {
		print(paste0("Flagging outliers for a dataset collected with slow sensor : ", dataset_id))
		rise_time = 0.05
		full_decay_time = 0.2
		
		if(toggle_short_peaks == TRUE){
					window = ceiling(minimum_window_duration/interFrame)
					enter = ceiling(0.75*rise_time/interFrame) 
					exit = ceiling(0.75*full_decay_time/interFrame) 
					
		} else {
		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2.5*rise_time/interFrame) 
		exit = ceiling(2.5*full_decay_time/interFrame) 
		}


	} else if (dataset_id == "GCaMP8f") {
		print(paste0("Flagging outliers for a dataset collected with medium-speed sensor : ", dataset_id))
		rise_time = 0.025
		full_decay_time = 0.150
		window = ceiling(minimum_window_duration/interFrame)
		enter = ceiling(2.5*rise_time/interFrame)
		exit = ceiling(2.5*full_decay_time/interFrame)

	}


if(dataset_id == "JF646") {
	print(paste0("Dataset ID is : ", dataset_id, ", we will take a percentile filter approach."))
		bin_groupers = c(groupers, "bins")
		data <- data_flagged %>%

		      	group_by_at(groupers) %>%
		      	mutate(bins = ntile(absoluteTime, 10)) %>%
		      	group_by_at(bin_groupers) %>%
		      	mutate(pctile = percent_rank(intensity_smooth),
					 alpha_str = case_when(flag == "Outlier" &  sensor == "JF646" ~ "isPeak",
									pctile > 0.3 & sensor == "JF646"  ~ "isPeak",
									pctile <= 0.3 & sensor == "JF646" ~ "NotPeak")) %>%
		      	ungroup() %>%
		      	group_by_at(groupers) %>%
		      	mutate(intensity_notPeak = ifelse(alpha_str == "NotPeak", intensity_smooth,NA),
						baseline_v0 = na.locf(intensity_notPeak,fromLast=FALSE,na.rm=FALSE),
						baseline=rollapply(baseline_v0,window,mean,fill=NA,align='center'),
						normIntensity = intensity_smooth/max(intensity_smooth,na.rm=TRUE),
		      			normBaseline = baseline/max(intensity_smooth,na.rm=TRUE),
		      			normStDev = signal_stdev/max(intensity_smooth,na.rm=TRUE),
		      				
		      			F = intensity_smooth,
		      			dF = (intensity_smooth - baseline),
		      			dFF = (intensity_smooth - baseline)/baseline,
		      				
		      			dFF_stdev = rollapply(data=dFF, width=window, FUN=sd,fill=NA, align="center"), 
		      			outlier_indices = case_when(flag == "Outlier" ~ "red")
		      			) 





} else {


data <- data_flagged %>%
      	group_by_at(groupers) %>%
      	mutate(normIntensity = intensity_smooth/max(intensity_smooth,na.rm=TRUE),
      				normBaseline = pseudoBaseline/max(intensity_smooth,na.rm=TRUE),
      				normStDev = signal_stdev/max(intensity_smooth,na.rm=TRUE),
      				
      				F = intensity_smooth,
      				dF = (intensity_smooth - pseudoBaseline),
      				dFF = (intensity_smooth - pseudoBaseline)/pseudoBaseline,
      				
      				dFF_idx = as.numeric(ifelse(flag=="notOutlier", dFF, "NA")),
    				dFF_signal_v0 = na.locf(dFF_idx, fromLast = FALSE, na.rm = FALSE),				
      				dFF_stdev = rollapply(data=dFF, width=window, FUN=sd,fill=NA, align="center"), 
      				outlier_indices = case_when(flag == "Outlier" ~ "red")
      				) 

}


drops <- c('dFF_idx',"dFF_signal_v0")  #extraneous data.frame columns clogging up memory
data<- data[ , !(names(data) %in% drops)]


print("Finished calculating dFF.")

toc()



tic()

source(paste0(path,"sub-functions/schmittTrig_v2.R"))
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/schmittTrig_v2.R") ## import Schmitt Trigger function with plotting


cumplus <- function(y) Reduce(function(a,b) a + b > 0, y, 0, accum=TRUE)[-1]

# threshold over which to identify a peak. Lower trigger = 1.5*sigma, upper trigger = 3.5*sigma
threshold = c(1.5,3.5)    
thresh_stdv = 1.5

peaksHeavy<-data %>%
		group_by_at(groupers) %>%
		do(schmittTrig(dataframe=., 
						time=.$absoluteTime, 
						dFF=.$dFF, 
						normIntensity=.$normIntensity, 
						normBaseline=.$normBaseline,
						#normSmooth = .$normSmooth, 
						normStDev=.$normStDev, 
						threshold=threshold, 
						std=.$dFF_stdev, 
						thresh_stdev=thresh_stdv)
					) %>% 	  
			mutate(#interFrame = absoluteTime - lag(absoluteTime),
					LongEntry = case_when(signal == 1 ~ 0,
											lead(signal, n= enter) == 1 ~ 1,
											signal == 0 ~ 0),
					LongExit = case_when(signal == 1 ~ 0,
											lag(signal, n= exit) == 1 ~ 1,
											signal == 0 ~ 0),
					#LongExitTest = 
					cumplus = cumplus(LongEntry - LongExit),
					temp = cumplus - c(0,pmin(0,diff(cumplus))),
					peak_idx = case_when(temp == 1 ~ "isPeak",
										temp == 0 ~ "notPeak")) %>%
			select(-temp) #get rid of extraneous columns

			

										 

print("Finished generously identifying likely peak indices using Schmitt Trigger.")
toc()




### UNCOMMENT THIS BLOCK OF CODE TO SAVE PLOTS FROM peakFinder_Heavy.R ######



# if(savePlots == TRUE) {

# 	library(ggeasy)
# 	library(scales)

# 	if(is.null(peaksHeavy$protocol)) {
# 		subDir <- paste0("peakFinderHeavyOutputs")	
# 	} else {
# 		subDir <- paste0("peakFinderHeavyOutputs",first(unique(peaksHeavy$protocol) ) )
  	
# 	}
	


# 	#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme_peakFinder_outputs.R") 
	

# 	#do things
# 	#want to save out schmitt Trigg 3 essentially

# 	.ROI_keyvars = rlang::syms(ROI_keys)
# 	tmp_df<- peaksHeavy %>% mutate(ROI_keys = paste(!!!.ROI_keyvars, sep="-"))

		
# 	if(exists('randomROIs')) {
# 		print("using the previous group of ROIs from peakFinder_Heavy.")
# 	} else {
# 		ROIs = unique(tmp_df$ROI_keys)
# 		ROIs_to_sample = ceiling( length(ROIs) /n ) 
# 		randomROIs<- sample(ROIs,ROIs_to_sample)

# 	}
	
# 	print( paste0("There are ", length(ROIs), " total unique ROIs in this dataset. We are going to save out traces for ", (1/n)*100, "% of them.") )

# 	 for (i in 1:length(randomROIs)) {

# 		tmp_data <- tmp_df %>%
#         dplyr::filter(ROI_keys == randomROIs[i]) 
#         lvl = mean(tmp_data$dFF_stdev, na.rm = TRUE)*threshold
#     	ymax = mean(tmp_data$max_dFF, na.rm=TRUE)
#     	ymin = mean(tmp_data$min_dFF, na.rm=TRUE)
   	






   
# 	## p1 - raw trace wwith outliers
# 	## p2 - normalized intensity with baseline overlay
# 	## p3 - dFF calculated with schmitTrig levels as geom_hline
# 	## p4 - signal binary

# 	p1<-ggplot(tmp_data, aes(x=absoluteTime, y = normIntensity)) +
#             geom_line(aes(fill='Normalized Signal'),size = 3) +
#             geom_point(aes(fill="Outliers Detected",colour=outlier_indices), size=5,shape=21,fill="white",stroke=2,alpha=0.5)+
#             labs( x="Time (s)",
#                   y="Normalized Intensity",
#                   title=unique(tmp_data$ROI_keys),
#                   subtitle="Normalized Trace with Outliers Identified")+
#             theme_tufte()+
#             scale_colour_identity(guide="legend")+
#             scale_y_continuous(labels = number_format(accuracy = 0.1))+
#             #coord_cartesian(ylim=c(ymin,ymax))+
#             #scale_colour_manual()+
#             scale_fill_manual(values=c('Normalized Signal'='black','Outliers Detected' = NA))+
#             my.theme+
#             theme(plot.title = element_text(size = 15),
#             		legend.position="none")+
#             easy_remove_x_axis()


#     p2<-ggplot(tmp_data, aes(x=absoluteTime, y = normIntensity)) +
#             geom_line(aes(colour="Normalized Signal"),size = 3 ) +
#             geom_line(aes(x=absoluteTime,y=normBaseline,colour="Estimated Baseline"),size=3.5) + 
            
#             labs( x="Time (s)",
#                   y="Normalized Intensity",
#                   title="",
#                   subtitle="Normalized Trace with Estimated Baseline")+
#             theme_tufte()+
#             scale_y_continuous(labels = number_format(accuracy = 0.1))+
            
#             scale_colour_manual(values=c('Normalized Signal'='black','Estimated Baseline'='red'))+
#             my.theme+
#             theme(legend.position="none")+
#             easy_remove_x_axis()   


#     p3<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF)) +
#             geom_line(aes(colour="Normalized Signal"),size = 3 ) +
#             geom_hline(colour="blue",yintercept=lvl,size=2) + 
            
#             labs( x="Time (s)",
#                   y=expression(Delta*"F/F"),
#                   title="",
#                   subtitle=expression(Delta*"F/F with Schmitt Trigger"))+
#             theme_tufte()+
#             coord_cartesian(ylim=c(ymin,ymax))+
#             scale_y_continuous(labels = number_format(accuracy = 0.1))+
#             scale_colour_manual(values=c('Normalized Signal'='black',"Schmitt Trigger" ='blue'))+
#             my.theme+
#             theme(legend.position="none")+
#             easy_remove_x_axis()    

#      p4<-ggplot(tmp_data, aes(x=absoluteTime, y = signal)) +
#              geom_line(aes(colour="Extracted Signal Indices"),size = 3 ) +
            
#              labs( x="Time (s)",
#                    y="Signal",
#                    title="",
#                    subtitle="")+
#              theme_tufte()+
#              scale_y_continuous(breaks=c(0.0,1.0),labels = number_format(accuracy = 0.1))+
#              coord_cartesian(ylim=c(-0.1,1.1))+

            
#              scale_colour_manual(values=c('Extracted Signal Indices'='red'))+
#              my.theme+
#              theme(legend.position="none",
#              		plot.title=element_blank())

#      fullplot=grid.arrange(p1, p2, p3,p4, nrow=4)



# 	ggsave(filename=paste0("peakFinderHeavy_trace#",i,".png"),plot=fullplot, device="png",dpi=600, units="in",width=12,height=20)
   



# rm(tmp_data)
# }
      

# }
# rm(tmp_df)


drops <- c("temp", "normIntensity", "normBaseline", "pseudoBaseline", "F", "dF", "dFF_stdev",'LongEntry','LongExit','cumplus')  #extraneous data.frame columns clogging up memory
peaksHeavy<- peaksHeavy[ , !(names(peaksHeavy) %in% drops)]


#warnings back on
options(warn = oldw)


rm(list=ls()[! ls() %in% c("peaksHeavy","groupers",'showGraphs','savePlots','mainDir','path', 'n','randomROIs','ROIs','ROIs_to_sample')])  #clear all vars in memory except for flagged data.

if(!is.null(dev.list())) dev.off()

