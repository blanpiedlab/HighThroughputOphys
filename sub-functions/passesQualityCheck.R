#passesQualityCheck.R
minTime = min(fitPeaks$absoluteTime, na.rm=TRUE)
maxTime = max(fitPeaks$absoluteTime, na.rm=TRUE)

initial_clean<-  suppressMessages(
                fitPeaks %>% group_by_at(groupers) %>% 
                #dplyr::filter(sum(signal) >= 3) %>%
                dplyr::filter(min(absoluteTime,na.rm=TRUE) != minTime & max(absoluteTime,na.rm=TRUE) != maxTime)
                )

getFitQuality<- suppressMessages( 
                initial_clean %>% group_by_at(groupers) %>%
                slice(which.max(dFF):length(dFF)) %>%
                summarise(fitQualityPass = all(.fitted == cummin(.fitted))) %>%
                dplyr::filter(fitQualityPass == TRUE)
                )

getAmplitude<- suppressMessages(
                initial_clean %>% group_by_at(groupers) %>% 
                summarise(amplitude = max(dFF)) %>% 
                ungroup() %>% 
                dplyr::filter(amplitude>= quantile(amplitude,0.025), amplitude <= quantile(amplitude, 0.975))
                )

getDuration<-   suppressMessages(
                initial_clean %>% group_by_at(groupers) %>% 
                summarise(duration = length(signal[which(signal == 1)])) %>% 
                ungroup() %>% 
                dplyr::filter(duration <= quantile(duration, 0.975))
                )

cleanPeaks <- suppressMessages(
                left_join(getAmplitude,getDuration, by = groupers)
                ) 
cleanPeaks <- suppressMessages(
                left_join(cleanPeaks,getFitQuality, by = groupers)
                )





fitsLabeled <- suppressMessages(
                left_join(fitPeaks, cleanPeaks)
                )   

rm(cleanPeaks,getDuration,getAmplitude,getFitQuality,initial_clean)



