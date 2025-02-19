#expFitter.R 
# get exponential decay fits for every peak identified in the dataset. 

tic()
library(ggforce)
source(paste0(path,"peakAnalysis/expFit.R")) ## import exponential fit function


fitData<-expDecayFit(df = df, groupers= groupers) #%>%
          #na.omit()

fitPeaks<- suppressMessages(left_join(df,fitData) %>%
                mutate(peakColor = ifelse(peakID == "NotPeak","noPeak", "isPeak"))
                )

fitPeaks$peakColor <- factor(fitPeaks$peakColor,
                                  levels = c("noPeak", "isPeak"))



rm(fitData, df_mutated)

print("Finished peak fitting routine.")
toc()