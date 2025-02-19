#plotRange.R
#Set variables for peak display plots using output from quality check
#set up variables for peakPlotter

 
peaksThatPass<- fitsLabeled %>% dplyr::filter(  fitQualityPass == TRUE, 
                                                    !is.na(amplitude), 
                                                    !is.na(duration)
                                                    )

windowed_ROIpeaks_ThatPass<- unique(peaksThatPass$windowed_ROIpeaks)
peaksToPlot<- windowed_ROIpeaks[windowed_ROIpeaks %in% windowed_ROIpeaks_ThatPass]


ymin= min(fitsLabeled$dFF, na.rm=TRUE)
ymax = max(peaksThatPass$amplitude, na.rm=TRUE)*1.1 



vids_sample<-as.numeric(length(vids))
ROIs_sample<- as.numeric(length(ROIs))
dirty_peakSample<- as.numeric(length(windowed_ROIpeaks))
cleaned_peakSample<- as.numeric(length(peaksToPlot))
