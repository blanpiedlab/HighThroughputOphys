#expFitter.R 
# get exponential decay fits for every peak identified in the dataset. 

tic()
library(ggforce)
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/expFit.R") ## import exponential fit function


fitData<-expDecayFit(df = data, groupers= groupers) #%>%
          #na.omit()

fitPeaks<- suppressMessages(left_join(data,fitData) %>%
                mutate(peakColor = ifelse(peakID == "NotPeak","noPeak", "isPeak"))
                )

fitPeaks$peakColor <- factor(fitPeaks$peakColor,
                                  levels = c("noPeak", "isPeak"))



rm(fitData)  #clear all vars in memory except for flagged data.

print("Finished peak fitting routine.")
toc()