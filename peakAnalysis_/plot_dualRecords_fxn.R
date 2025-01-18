#plot_dualRecord_minis.R

#Written by Samuel T. Barlow 
# 11.8.22



######## for every glutamate event, is there a following calcium event? ########## 




#Plotting dataset from combined imaging

library(ggforce)

#subDir <- "dualRecordings_v7"

#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 

library(scales)



### This function is a function-ification of an existing script. Since I will be embedding it in an existing script, it should accept existing inputs. 

### plot_dualRecords should accept a cleanedup list of traces, perform the necessary dataframe manipulations, and subset them accordingly. 
### for now, save all the data 



plot_dualRecords<- function( subsetted_traces,interspike_threshold,tag){


minpositive = function(x) min(x[x > 0 ])

        tmp_trace <- subsetted_traces 
        graphTitle = first(unique(tmp_trace$trackROI_key) )



                    slicePeaks <- tmp_trace %>% 
                                      dplyr::filter(peakID != "NotPeak") %>% 
                                      group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
                                      slice(which.max(dFF)) %>% 
                                      ungroup() %>% 
                                      group_by(chemical_condition,trackROI_key) %>% 
                                      arrange(absoluteTime) %>%
                                      mutate(
                                             interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
                                             isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
                                             timeCheck = ifelse(isValid == TRUE & interSpike < interspike_threshold, TRUE, FALSE),
                                              isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
                                                                          sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
                                                                          sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
                                              checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)
                                               )

                    checkInterspike<- slicePeaks %>% dplyr::filter(timeCheck == TRUE)
                    getOrphans<- slicePeaks %>% dplyr::filter(timeCheck == FALSE)


                    sliceGluPeaks_v2<- tmp_trace %>%  dplyr::filter(peakID != "NotPeak") %>% 
                                                          group_by(sensor,trackROI_key,chemical_condition,peakID) %>% 
                                                          slice(which.max(dFF)) %>% 
                                                          ungroup() %>% 
                                                          group_by(chemical_condition,trackROI_key) %>% 
                                                          arrange(absoluteTime) %>%
                                                          mutate(interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
                                                                    isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
                                                                     timeCheck = ifelse(isValid == TRUE & interSpike < interspike_threshold, TRUE, FALSE),
                                                                      isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
                                                                                                  sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
                                                                                                  sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
                                                                      checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE)   ) %>%
                                                        dplyr::filter(sensor == "GluSnFR3")

                    get_Glu_peakPositions = sliceGluPeaks_v2                                      
                    sliceGluPeaks_v2$timeCheck[is.na(sliceGluPeaks_v2$timeCheck)] <- FALSE
                    sliceGluPeaks_v2$sensor[sliceGluPeaks_v2$sensor == "GluSnFR3"] <- "JF646"
                    sliceGluPeaks_v2 <- sliceGluPeaks_v2 %>% ungroup() %>% group_by(chemical_condition) %>%mutate(yposX = 1.0*ifelse(isTransmission == FALSE, 0.9, NA),
                                                                                yposO = 1.0*ifelse(isTransmission == TRUE, 0.9, NA))
                    sliceTime <- sliceGluPeaks_v2 %>% select(absoluteTime) %>% mutate(sensor = "JF646")





facet.labels <- c("iGluSnFR3","JF646-BAPTA-AM")
names(facet.labels)<- c("GluSnFR3", "JF646")
tmp_trace$chemical_condition = factor(tmp_trace$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Vehicle", "AP5"))
sliceTime$chemical_condition = factor(sliceTime$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Vehicle", "AP5"))
sliceGluPeaks_v2$chemical_condition = factor(sliceGluPeaks_v2$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Vehicle", "AP5"))
get_Glu_peakPositions$chemical_condition = factor(get_Glu_peakPositions$chemical_condition,levels=c("ctrl","vehicle","APV"), labels=c("Pre-Treatment", "Vehicle", "AP5"))


 
 custom_y <- tribble(
  ~sensor, ~absoluteTime, ~dFF,
  "GluSnFR3", 0, -0.2,
  "GluSnFR3", 0, 5,
  "JF646", 0, -0.2,
  "JF646", 0, 1.0
                    )

 my_y_breaks <- function(x) {if(max(x) <= 1.1)  seq(0.0,1.0,by=0.5) else if(max(x) <= 5 & max(x) >= 1.2) seq(0.0,5,by=2.5) else seq(0.0,7.5,by=2.5)}
 my_breaks <- function(x) { if(max(x) < 20)  seq(0,15,5) else seq(0,25,5) }

 tmp_trace$sensor = as.factor(as.character(tmp_trace$sensor))
 levels(tmp_trace$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-AM"))
 
 sliceTime$sensor = as.factor(as.character(sliceTime$sensor))
 levels(sliceTime$sensor) <- c(expression(JF[646]*"-BAPTA-AM"))

 sliceGluPeaks_v2$sensor = as.factor(as.character(sliceGluPeaks_v2$sensor))
 levels(sliceGluPeaks_v2$sensor) <- c(expression(JF[646]*"-BAPTA-AM"))

 get_Glu_peakPositions$sensor = as.factor(as.character(get_Glu_peakPositions$sensor))
 levels(get_Glu_peakPositions$sensor) <- c("iGluSnFR3")

 custom_y$sensor = as.factor(as.character(custom_y$sensor))
 levels(custom_y$sensor) <- c("iGluSnFR3",expression(JF[646]*"-BAPTA-AM"))


tracePlot<-ggplot(tmp_trace, aes(x=absoluteTime, y = dFF, colour = sensor, group = 1)) +
            geom_vline(data=sliceTime, aes(xintercept = absoluteTime), colour = "black", alpha = 0.5, linetype = "longdash", size = 2) +
            geom_text(data=sliceGluPeaks_v2,aes(label = "X",x = absoluteTime, y = yposX), colour="blue", size= 14,alpha=0.75)+
            geom_text(data=sliceGluPeaks_v2,aes(label = "O",x = absoluteTime, y = yposO), colour="black", size= 14,alpha=0.75)+
            geom_blank(data=custom_y,aes(absoluteTime,dFF))+
            geom_line(size = 2) + ####the actual data
            geom_point(data=get_Glu_peakPositions, aes(x=absoluteTime,y=dFF), shape=23, fill="blue",alpha=0.7,colour="blue", size=8)+
            
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  #title="",
                  #subtitle=graphTitle,
                  tag = tag)+
            scale_colour_manual(name = "", values = c("forestgreen", "#FF33FF"))+
            theme_tufte()+
            my.theme+
             theme(legend.position="none",strip.text.y=element_blank())+
            facet_grid(sensor~chemical_condition, labeller = labeller(sensor = facet.labels),scales="free_y")+
            scale_y_continuous(labels = number_format(accuracy = 0.01), breaks=my_y_breaks, expand = c(0.15,0))+
            scale_x_continuous(expand=c(0,0),breaks = my_breaks)#+
           # coord_cartesian(ylim=c(-0.1,1.5))
             


ggsave(filename=paste0(graphTitle,"dualRecording.png"),plot=tracePlot, bg="white",device="png",dpi=600, units="in",width=30,height=20)

tracePlot    
}