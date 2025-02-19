#getPeakLabels.R
colnameGetter = function(d,...){
  #browser()
  vars  <- quos(...)
  d %>% mutate(ROIpeak = paste(!!!vars, sep="_"))
}

label0= "ROINumber"
label1= "peakID"
label2= "windowedPeakID"

forVid = c(label0,label1,label2)
vid_groupers = setdiff(groupers,forVid)
getVidLabels<- peaks %>% ungroup() %>% select(one_of(vid_groupers)) 
getVidLabels <- unique(getVidLabels)

vids <- apply(getVidLabels[, vid_groupers], 1, paste, collapse = '_')



forROI = c(label1,label2)
ROI_groupers = setdiff(groupers,forROI)
getROILabels<- peaks %>% ungroup() %>% select(one_of(ROI_groupers)) 
getROILabels <- unique(getROILabels)

ROIs <- apply(getROILabels[, ROI_groupers], 1, paste, collapse = '_')

ROIpeak_groupers = setdiff(groupers,label2)
getPeakLabels <- peaks %>% dplyr::filter(peakID != "NotPeak") %>% select(one_of(ROIpeak_groupers)) 
getPeakLabels <- unique(getPeakLabels)

ROIpeaks <- apply(getPeakLabels[, ROIpeak_groupers], 1, paste, collapse = '_')


windowed_peakGroupers = setdiff(groupers,label1)
getWindowed_PeakLabels <- peaks %>% dplyr::filter(windowedPeakID != "NotPeak") %>%  select(one_of(windowed_peakGroupers))  
getWindowed_PeakLabels <- unique(getWindowed_PeakLabels)

windowed_ROIpeaks <- apply(getWindowed_PeakLabels[, windowed_peakGroupers], 1, paste, collapse = '_')




peaksLabeled<- peaks %>% group_by_at(groupers) %>% unite("windowed_ROIpeaks",all_of(windowed_peakGroupers), remove=FALSE) %>% unite("ROIpeaks",all_of(ROIpeak_groupers), remove=FALSE) %>% unite("ROIs", all_of(ROI_groupers), remove=FALSE)


rm(list=ls()[! ls() %in% c("peaks","peaksCleaned","peaksLabeled","groupers","ROIpeaks","ROIs","windowed_ROIpeaks","vids")])  #clear all vars in memory except for flagged data.

