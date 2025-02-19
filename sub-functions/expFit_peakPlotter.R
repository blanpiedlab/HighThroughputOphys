
#suppress warnings
oldw <- getOption("warn")
options(warn = -1)





#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R

subDir <- "peakFitters"   


source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/def_stim_vlines.R")


library(ggforce)



#note: paradigm will break when peaks are spaced by less than 50 frames
# paradigm breaks because sampling the dplyr::filter call will only take the frames associated with the peak. i.e. PeakID possess only 2 frames, only 2 frames long graph.
# let's just test.  






randomROIpeaks<- sample(peaksToPlot,n)
interFrame = as.vector(fitsLabeled %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)


tic()
#for (i in 1:length(ROIs)) {
  
  # Establish a tmp_data vector that can be used plot only the ROIs from a specific region for each iteration of the outer loop. 
  
   
  for (i in 1:length(randomROIpeaks)) {
  suppressMessages(tmp_data <- left_join(sync_async_fitsLabeled,halfWidths) %>%
        dplyr::filter(windowed_ROIpeaks == randomROIpeaks[i]) %>%
        mutate(timeRef = paste(!!!.keyvars, sep="-")  )
        )

    tmp_ID = randomROIpeaks[i]
    tmp_stim_ID = unique(tmp_data$timeRef)  
    stimParadigm_times = def_stim_vlines(tmp_stim_ID)

    vline_switch = !is.null(stimParadigm_times)

    xmin = min(tmp_data$absoluteTime)-(10*interFrame)
    xmax = max(tmp_data$absoluteTime)+(10*interFrame)
    #ymax = max(fitPeaks$dFF)
    #ymin = min(fitPeaks$dFF)*0.5
    #for (j in 1:length(peakPage)) {
        tmp_tracePlot<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF, colour = factor(peakColor), group =1)) +
            geom_line(size = 1.5) + 
            geom_path(aes(x=absoluteTime, y = .fitted, group=1,linetype="expFit"), colour="red",size=2,alpha=0.7)+
            geom_point(aes(x=riseEdgeTime,y=riseEdgeAmplitude), colour = "green", size=4)+                                                                                      #where is riseEdge
            geom_segment(aes(x = rise_half_time, y = rise_half_amplitude, xend = decay_half_time, yend = decay_half_amplitude, ),colour = "blue",lty="solid",size=1.5)+        #define segment showing half-width
            geom_point(aes(x=rise_half_time,y=rise_half_amplitude),colour="blue",size=4)+                                                                                       #where is rise_half_point
            geom_point(aes(x=decay_half_time,y=decay_half_amplitude),colour="blue",size=4)+                                                                                     #where is decay_half_point
            geom_point(aes(x=maxTime,y=maxAmplitude),colour="orange",size=4)+                                                                                                   #where is peak maximum
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  title="",
                  subtitle=ifelse(!is.null(stimParadigm_times), paste0("First stim estimated at:",stimParadigm_times[1], " seconds"), NA))+
            scale_colour_manual(values = c("noPeak" = "lightgrey", "isPeak" ="black"))+
            scale_linetype_manual(name="Exp Decay Fit", values="solid")+
            theme_tufte()+
            coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))+
            my.theme+
            theme(legend.position = "right",
                  legend.justification = c('right','top'),
                  #legend.position = "none",
                  #legend.title = element_text(colour="black", size=20, family="Serif"), 
                  #legend.text=element_text(colour="black", size=16, family="Serif"),
                  axis.text.x=element_text(colour="black", size=20, family="Serif"),
                  axis.text.y=element_text(colour="black", size=20, family="Serif"),
                  axis.title=element_text(colour="black", size=24, family="Serif"),
                  axis.ticks.length=unit(.25, "cm"),
                  axis.ticks = element_line(size=1),
                  strip.text = element_text(colour="black", size = 16, family = "Serif"))+
            
            facet_wrap(~windowed_ROIpeaks, scales ="free",ncol=1,nrow=1)
      
    
   ggsave(filename=paste0("peakFits_",i,".png"),plot=tmp_tracePlot, device="png",dpi=300, units="in",width=12,height=12)
    

   tmp_tracePlot<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF, colour = factor(timeClass), group =1)) +
            geom_line(size = 1.5) + 
            {if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "black", alpha = 0.4, linetype = "longdash", size = 1)} +
                        

            #geom_path(aes(x=absoluteTime, y = .fitted, group=1,linetype="expFit"), colour="red",size=2,alpha=0.7)+
            #geom_point(aes(x=riseEdgeTime,y=riseEdgeAmplitude), colour = "green", size=4)+                                                                                      #where is riseEdge
            #geom_segment(aes(x = rise_half_time, y = rise_half_amplitude, xend = decay_half_time, yend = decay_half_amplitude, ),colour = "blue",lty="solid",size=1.5)+        #define segment showing half-width
            #geom_point(aes(x=rise_half_time,y=rise_half_amplitude),colour="blue",size=4)+                                                                                       #where is rise_half_point
            #geom_point(aes(x=decay_half_time,y=decay_half_amplitude),colour="blue",size=4)+                                                                                     #where is decay_half_point
            #geom_point(aes(x=maxTime,y=maxAmplitude),colour="orange",size=4)+                                                                                                   #where is peak maximum
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  title="",
                  subtitle=ifelse(!is.null(stimParadigm_times), paste0("First stim estimated at:",stimParadigm_times[1], " seconds"), NA))+
            scale_colour_manual(values = c("noPeak" = "lightgrey", "sync" ="blue","async" = "red", "spontaneous" = "orange"))+
            scale_linetype_manual(name="Exp Decay Fit", values="solid")+
            theme_tufte()+
            coord_cartesian(xlim=c(xmin,xmax), ylim=c(ymin,ymax))+
            my.theme+
            theme(legend.position = "right",
                  legend.justification = c('right','top'),
                  #legend.position = "none",
                  #legend.title = element_text(colour="black", size=20, family="Serif"), 
                  #legend.text=element_text(colour="black", size=16, family="Serif"),
                  axis.text.x=element_text(colour="black", size=20, family="Serif"),
                  axis.text.y=element_text(colour="black", size=20, family="Serif"),
                  axis.title=element_text(colour="black", size=24, family="Serif"),
                  axis.ticks.length=unit(.25, "cm"),
                  axis.ticks = element_line(size=1),
                  strip.text = element_text(colour="black", size = 16, family = "Serif"))+
            
            facet_wrap(~windowed_ROIpeaks, scales ="free",ncol=1,nrow=1)
      
    
   ggsave(filename=paste0("peakFits_",i,"_sync_async_classified.png"),plot=tmp_tracePlot, device="png",dpi=300, units="in",width=12,height=12)
    











    print(paste0("Finished saving classified peaks for: ", randomROIpeaks[i]))
  #}    
  
  rm(list=ls(pattern="tmp_"))
}

toc()

rm(i,n,"peaksToPlot","windowed_ROIpeaks_ThatPass")


#warnings back on
options(warn = oldw)

# Post-script: This script is intended to guide the analysis of laser-power vs. stimulation number for improving our understanding of how GluSnFR behaves under different conditions. 


