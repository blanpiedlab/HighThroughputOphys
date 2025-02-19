#suppress warnings
oldw <- getOption("warn")
options(warn = -1)



subDir <- "traceFitters"  


#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/def_stim_vlines.R")

library(ggforce)


  





randomROIs<- sample(ROIs,n)

tic()
#for (i in 1:length(ROIs)) {
  
  # Establish a tmp_data vector that can be used plot only the ROIs from a specific region for each iteration of the outer loop. 
  .keyvars = rlang::syms(keys)
   
  for (i in 1:length(randomROIs)) {
  tmp_data <- sync_async_fitsLabeled %>%
        dplyr::filter(ROIs == randomROIs[i]) %>%
        mutate(timeRef = paste(!!!.keyvars, sep="-") )
    tmp_ID = randomROIs[i]


    tmp_stim_ID = unique(tmp_data$timeRef)  
    stimParadigm_times = def_stim_vlines(tmp_stim_ID)

    vline_switch = !is.null(stimParadigm_times)


    #ymax = max(fitsCleaned$dFF)
    #ymin = min(fitsCleaned$dFF)

        tmp_tracePlot<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF, colour = factor(peakColor), group =1)) +
            geom_line(size = 1) + 
            geom_path(aes(x=absoluteTime, y = .fitted, group=1,linetype="expFit"), colour="red",size=2,alpha=0.7)+
            {if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "black", alpha = 0.4, linetype = "longdash", size = 1)} +
            
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  title="",
                  subtitle=ifelse(!is.null(stimParadigm_times), paste0("First stim estimated at:", stimParadigm_times[1], " seconds"), NA))+
            scale_colour_manual(values = c("noPeak" = "lightgrey", "isPeak" ="black"))+
            scale_linetype_manual(name="Exp Decay Fit", values="solid")+
            theme_tufte()+
            coord_cartesian(ylim=c(ymin,ymax))+
            my.theme+
            theme(legend.position = "right",
                  legend.justification = c('right','top'),
                  #legend.position = "none",
                  #legend.title = element_text(colour="black", size=20, family="Serif"), 
                  #legend.text=element_text(colour="black", size=16, family="Serif"),
                  axis.text.x=element_text(colour="black", size=30, family="Serif"),
                  axis.text.y=element_text(colour="black", size=30, family="Serif"),
                  axis.title=element_text(colour="black", size=36, family="Serif"),
                  axis.ticks.length=unit(.25, "cm"),
                  axis.ticks = element_line(size=1),
                  strip.text = element_text(colour="black", size = 24, family = "Serif"))+
            
            facet_wrap(~ROIs, scales ="free",ncol=1,nrow=1)
      
    
   ggsave(filename=paste0("traceFits_",i,".png"),plot=tmp_tracePlot, device="png",dpi=300, units="in",width=30,height=15)

   ##testing sync/async classification
   tmp_tracePlot<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF, colour = factor(timeClass), group =1)) +
            geom_line(size = 1) + 
            #geom_path(aes(x=absoluteTime, y = .fitted, group=1,linetype="expFit"), colour="black",lty="dashed",size=1,alpha=0.7)+
            {if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "black", alpha = 0.4, linetype = "longdash", size = 1)} +
            
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  title="",
                  subtitle=ifelse(!is.null(stimParadigm_times), paste0("First stim estimated at:",stimParadigm_times[1], " seconds"), NA))+
            scale_colour_manual(values = c("noPeak" = "lightgrey", "sync" ="blue","async" = "red", "spontaneous" = "orange"))+
            scale_linetype_manual(name="Exp Decay Fit", values="solid")+
            theme_tufte()+
            coord_cartesian(ylim=c(ymin,ymax))+
            my.theme+
            theme(legend.position = "right",
                  legend.justification = c('right','top'),
                  #legend.position = "none",
                  #legend.title = element_text(colour="black", size=20, family="Serif"), 
                  #legend.text=element_text(colour="black", size=16, family="Serif"),
                  axis.text.x=element_text(colour="black", size=30, family="Serif"),
                  axis.text.y=element_text(colour="black", size=30, family="Serif"),
                  axis.title=element_text(colour="black", size=36, family="Serif"),
                  axis.ticks.length=unit(.25, "cm"),
                  axis.ticks = element_line(size=1),
                  strip.text = element_text(colour="black", size = 24, family = "Serif"))+
            
            facet_wrap(~ROIs, scales ="free",ncol=1,nrow=1)
      
    
   ggsave(filename=paste0("traceFits_",i,"sync_async_classified.png"),plot=tmp_tracePlot, device="png",dpi=300, units="in",width=30,height=15)
    

   if(vline_switch == FALSE) {
     print("Trace belongs to Sham stimulus group.")    
     print(paste0("Finished saving classified traces for: ", randomROIs[i]))
     rm(list=ls(pattern="tmp_"))
     next
    }
   

    if(vline_switch == TRUE & stimParadigm_times[1] == 0.45){

    print("The stimList for this dataset used default.onset.")    
    print(paste0("Finished saving classified traces for: ", randomROIs[i]))
    rm(list=ls(pattern="tmp_"))
    next
    }    
    
    if(vline_switch == TRUE & stimParadigm_times[1] != 0.45){

    print("The stimList for this dataset used a derived find.onset.")    
    print(paste0("Finished saving classified traces for: ", randomROIs[i]))
    rm(list=ls(pattern="tmp_"))
    next
    }
   



    rm(list=ls(pattern="tmp_"))
}

toc()

rm(i,n)



#warnings back on
options(warn = oldw)

# Post-script: This script is intended to guide the analysis of laser-power vs. stimulation number for improving our understanding of how GluSnFR behaves under different conditions. 


