oldw <- getOption("warn")
options(warn = -1)



subDir <- paste0("peakFitters_v1" ) 
  


#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/def_stim_vlines.R")

library(ggforce)

#tracePlotter is nice, in that it doesn't require too much to run. Here we can take bits and pieces from getPeakLabels.R to get it to run. 

if (hide_titles == TRUE) {
    title_switch = TRUE
} else {
    title_switch=FALSE
}


ROIs = unique(df$ROI_key)
ROIs_to_sample = ceiling( length(ROIs) /n ) 
randomROIs<- sample(ROIs,ROIs_to_sample)

print( paste0("There are ", length(ROIs), " total unique ROIs in this dataset. We are going to save out traces for ", (1/n)*100, "% of them.") )

tic()
  
  # Establish a tmp_data vector that can be used plot only the ROIs from a specific region for each iteration of the outer loop. 
  #.keyvars = rlang::syms(keys)
scalefactor = 0.75   
trace_width = scalefactor*16
trace_height=scalefactor*16


  for (i in 1:length(randomROIs)) {
        tmp_data <- df %>%
            dplyr::filter(ROI_key == randomROIs[i]) 
        
        newPeakColors = list(Signal = "isPeak", Baseline = "noPeak")
        levels(tmp_data$peakColor) <- newPeakColors
        print(paste0("ROI_key is: ",unique(tmp_data$ROI_key) ) )
            
        tmp_peaks_total <- unique(tmp_data$windowedPeakID)
        print("Found these windowedPeakIDs:")
        print(tmp_peaks_total)

        for (j in 1:length(tmp_peaks_total)) {

            tmp_peak <- tmp_data %>%
                        dplyr::filter(windowedPeakID == tmp_peaks_total[j]) %>%
                        mutate(windowed_normTime = absoluteTime - min(absoluteTime))
            
            if(add_time_dots == TRUE){
                    tmp_time_dots <- output_data %>% 
                                    dplyr::filter(ROI_key == randomROIs[i], 
                                                    windowedPeakID == tmp_peaks_total[j])
    
            }

            
            tmp_ymax = max(tmp_peak$dFF)
            if(is.finite(tmp_ymax) & tmp_ymax < 0.3) {
                next
            }

            if(is.finite(tmp_ymax) & tmp_ymax > 1.5) {
                ymax = tmp_ymax
            } else {
                ymax = 1.5
            }
            ymin = -0.1
            xmin = min(tmp_peak$absoluteTime)
            xmax = xmin+0.5
            tmp_ID = paste0(randomROIs[i],"--",unique(tmp_peak$windowedPeakID))
            if(add_max == TRUE) {
                    #usual_frame_time<-median(tmp_data$interFrame,na.rm=TRUE)
                    #dur_TTL = 0.01 #s

                    
                    #TTL_index = unique(tmp_data$TTL_start)
                    max_pt<-tmp_data %>% slice(which.max(dFF)) %>% dplyr::filter(peakID != "NotPeak")
                    pct90_line <- as.numeric(0.9*max_pt$dFF)
                    pct10_line <- as.numeric(0.1*max_pt$dFF)
                    if(length(!is.na(max_pt$dFF)) > 0 ){ tmp_max = 1} else { tmp_max = 0}
                }

            #check variable number
            check_cols<- c('half_rise_time','half_rise_dFF',
                            'half_decay_time','half_decay_dFF',
                            'rise10_time','rise10_dFF',
                            'rise90_time','rise90_dFF',
                            'decay10_time','decay10_dFF',
                            'decay90_time','decay90_dFF')
            
            tmp_peakColors<- length(unique(tmp_peak$peakColor) )
            if('.fitted' %in% colnames(tmp_peak)){tmp_decay = 1} else {tmp_decay = 0}
            
            if(all(check_cols %in% colnames(tmp_time_dots)) ){tmp_time_vars = 3} else {tmp_time_vars = 0}

            sum_vars = sum(tmp_peakColors,tmp_decay,tmp_time_vars, tmp_max)

            if(sum_vars < 7) {
                print(paste0("Sum of the variables for legend is: ", sum_vars))
                check_tmp_times <- tmp_time_dots[,names(tmp_time_dots) %in% check_cols]
                print(colnames(check_tmp_times))    
                print("skipping to next loop")
                next
            }
            
            parse.labels <- function(x) parse(text = x)

            tmp_tracePlot<-ggplot(tmp_peak, aes(x=absoluteTime, y = dFF, colour = peakColor, group =1)) +
                
                            
                            # {if(add_scalebars)geom_segment(data=scalebar_x, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 3, colour="black")}+
                            # {if(add_scalebars)geom_segment(data=scalebar_y, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 3, colour="black")}+
                            # {if(add_scalebars)geom_label(aes(x=label_x$x,y=label_x$y, label=label_x_expr), nudge_y = (user_y*-0.15), size= 7, colour="black",parse=TRUE,label.size=NA,family="sans")}+
                            # {if(add_scalebars)geom_label(aes(x=label_y$x,y=label_y$y, label=label_y_expr), nudge_x = (scalebar_x_length*-0.7), size= 7, colour="black",parse=TRUE,label.size=NA,family="sans")}+
                            
                           
                            {if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=half_rise_time,y=half_rise_dFF,xend=half_decay_time,yend=half_decay_dFF, colour="t[1/2]"), size= 2, alpha=1)}+ #segment to join the two
                            geom_path(size = 1.5,lineend='round',alpha=1.0) + 
                            geom_path(aes(x=absoluteTime, y = .fitted, group=1, colour="Exp.~Decay~Fit"),linetype="solid", size=2,alpha=0.8)+
                            #{if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "grey21", alpha = 1, linetype = "longdash", size =0.9)} +
                            #{if(vline_switch)geom_label(aes(x = stimParadigm_times, y = ifelse(is.finite(ymax), ymax*1.24, 0 ) ),label = "Stimulus", colour = "black", family='sans', size=7,label.size=NA)}+
                            
                             #half_width
                            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = half_rise_time, y = half_rise_dFF, colour="t[1/2]"), shape=21, fill="aquamarine4",size= 6, alpha=0.85)}+ #half_t[rise]
                            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = half_decay_time, y = half_decay_dFF, colour="t[1/2]"), shape=21, fill="aquamarine4",size= 6, alpha=0.85)}+ #half_decay_time


                             ##dashed lines indicating 90 and 10% dFF levels
                           
                            #{if(add_time_dots)geom_segment(aes(x=max_pt$absoluteTime-(0.5*scalebar_x_length), xend=max_pt$absoluteTime+(scalebar_x_length), y=pct90_line,yend=pct90_line), colour = "grey41", alpha = 1, linetype = "dashed", size = 0.6)}+
                            #{if(add_time_dots)geom_segment(aes(x=max_pt$absoluteTime-(0.5*scalebar_x_length), xend=decay10_lineend+(scalebar_x_length), y=pct10_line,yend=pct10_line), colour = "grey41", alpha = 1, linetype = "dashed", size = 0.6)}+
                            #{if(add_time_dots)geom_label(aes(x = max_pt$absoluteTime+(scalebar_x_length), y = pct90_line), nudge_x = 0.035, label = "90%",fontface="italic", colour = "black", family='sans', size=6.5,label.size=NA,alpha=0.95,fill=NA)}+
                            #{if(add_time_dots)geom_label(aes(x = tmp_time_dots$decay10_time+(scalebar_x_length), y = pct10_line), nudge_x = 0.035,label = "10%",fontface="italic", colour = "black", family='sans', size=6.5,label.size=NA,alpha=0.95,fill=NA)}+

                            
                            # #rise_time
                             #{if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=rise10_time,y=rise10_dFF,xend=rise90_time,yend=rise90_dFF), size= 3, alpha=0.75,colour="green")}+ #segment to join the two
                             {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = rise10_time, y = rise10_dFF, colour="t[rise]"), shape=21, fill="deepskyblue",size= 5, alpha=0.85)}+ #10riseTime
                             {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = rise90_time, y = rise90_dFF, colour="t[rise]"), shape=21, fill="deepskyblue",size= 5, alpha=0.85)}+ #90riseTime

                            # #  {if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=decay10_time,y=decay10_dFF,xend=decay90_time,yend=decay90_dFF), size= 3, alpha=0.75,colour="orange")}+ #segment to join the two
                             {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = decay10_time, y = decay10_dFF, colour="t[decay]"), shape=21, fill="darkorchid", size= 5, alpha=0.85)}+ #10decayTime
                             {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = decay90_time, y = decay90_dFF,colour="t[decay]"), shape=21, fill="darkorchid",size= 5, alpha=0.85)}+ #90decayTime


                            #decay_time

                            
                           
                         
                            {if(add_max)geom_point(data=max_pt, aes(x=absoluteTime,y=dFF, colour = "Peak[Delta*F/F]"), shape=23, fill="blue",alpha=0.65, size=8)}+
                                                  
                           
                           
                            labs( x="Time (s)",
                                  y=expression(Delta*"F/F"),
                                  subtitle=tmp_ID)+
                            scale_colour_manual(breaks = c("Signal", "Baseline","Exp.~Decay~Fit", "t[rise]", "t[decay]", "t[1/2]", "Peak[Delta*F/F]"), 
                                                 values = c( "Signal" ="black","Baseline" = "grey72","Exp.~Decay~Fit"="red", "t[rise]" = "deepskyblue", "t[decay]" = "darkorchid","t[1/2]" = "aquamarine4", "Peak[Delta*F/F]" = "blue"),
                                                 labels = parse.labels)+

                            #scale_linetype_manual(breaks = c("Signal","Baseline",  "Exp Decay Fit"),
                            #                        values = c( "Signal" = "solid", "Baseline" = "solid","Exp Decay Fit" = "solid") )+
                            theme_tufte()+

                            guides(colour = guide_legend(override.aes = list(linetype = c(1,1,1,NA,NA,1,NA),
                                                                                colour = c("black","grey72","red", "deepskyblue","darkorchid","aquamarine4","blue"),
                                                                                shape = c(NA,NA,NA,21,21,21,23), 
                                                                                fill = c("black","grey72","red","deepskyblue","darkorchid","aquamarine4","blue"),
                                                                                alpha = c(1,1,0.9,0.85,0.85,0.85,0.65),
                                                                                size = c(2,2,2,5,5,6,8)
                                                                                )
                                                                            )
                                                                        )+
                            coord_cartesian(ylim=c(ymin,ymax),xlim=c(xmin,xmax))+
                            my.theme+
                            theme(plot.subtitle = element_text(size=16,family="sans"),
                                                                legend.position = c(0.8,0.8),
                                                                legend.text.align = 0,
                                                                legend.spacing.y = unit(0, "lines"),
                                                                legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                                    )
                            # {if(title_switch)theme(plot.title=element_text(size=14, family="sans"),
                            #                         legend.position = c(0.8,0.8),
                            #                         legend.text.align = 0,
                            #                         #legend.spacing.y = unit(0, "lines"),
                            #                         #legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                            #                         #)#element_blank(), 
                            #                             plot.subtitle=element_blank(),
                            #                             strip.text=element_blank(),
                            #                             plot.margin = margin(1,1,1,1, "cm"),
                            #                             axis.line = element_blank(),
                            #                             axis.text.x = element_blank(),
                            #                             axis.text.y = element_blank(),
                            #                             axis.title = element_blank(),
                            #                             axis.ticks = element_blank())}+
                           
                            #facet_wrap(~ROI_key,ncol=1,nrow=1)
                      
                    
                   ggsave(filename=paste0("peakFits_trace",i,"peak",j,"_.jpeg"),plot=tmp_tracePlot, device="jpeg",bg='white',dpi=300, units="in",width=trace_width,height=trace_height)

  
 


    rm(tmp_peak)
    }
       rm(list=ls(pattern="tmp_"))
}

toc()

rm(i,n,ROIs,randomROIs,ROIs_to_sample,plotChars,xmax,ymin,ymax,vline_switch,title_switch,trace_width,trace_height)



#warnings back on
options(warn = oldw)

