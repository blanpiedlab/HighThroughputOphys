#cleanupEnv.R
#run after peakPlotter

print("The total number of videos analyzed : ")
print(vids_sample)
print("The total number of ROI intensity-time traces analyzed : ")
print(ROIs_sample)
print("The total number of peaks identified : ")
print(dirty_peakSample)
print("The total number of peaks remaining after quality control check : ")
print(cleaned_peakSample)





rm(list=ls()[! ls() %in% c("fitsCleaned", "halfWidths", "mainDir",'groupers','groupers_plus','stimList')])  #clear all vars in memory except for flagged data.
