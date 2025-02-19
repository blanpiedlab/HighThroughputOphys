#suppress warnings
oldw <- getOption("warn")
options(warn = -1)



subDir <- paste0("traceFitters_v1Ca" ) 
  


#set figure directory with mainDir and subDirs articulated in run_pipeline_v3.R
source(paste0(path,"sub-functions/setFigureDirectory.R")) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"sub-functions/def_stim_vlines.R"))
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))

library(ggforce)

#tracePlotter is nice, in that it doesn't require too much to run. Here we can take bits and pieces from getPeakLabels.R to get it to run. 

if (hide_titles == TRUE) {
    title_switch = TRUE
} else {
    title_switch=FALSE
}

if(pub_quality ==TRUE) {
    dpi_val = 900
} else {
    dpi_Val = 300
}
special_check<-peak_stats %>% dplyr::filter(Ca == "1Ca")

ROIs = unique(special_check$ROI_key)
ROIs_to_sample = ceiling( length(ROIs) /n ) 
randomROIs<- sample(ROIs,ROIs_to_sample)

print( paste0("There are ", length(ROIs), " total unique ROIs in this dataset. We are going to save out traces for ", (1/n)*100, "% of them.") )

tic()
  
  # Establish a tmp_data vector that can be used plot only the ROIs from a specific region for each iteration of the outer loop. 
  #.keyvars = rlang::syms(keys)
   
  for (i in 1:length(randomROIs)) {
  tmp_data <- df %>%
        dplyr::filter(ROI_key == randomROIs[i]) 
  
  if(add_time_dots == TRUE){
  tmp_time_dots <- peak_stats %>% 
        dplyr::filter(ROI_key == randomROIs[i])
    
    }
   


        xmax = max(tmp_data$absoluteTime,na.rm=TRUE)
        ymax = max(tmp_data$dFF,na.rm=TRUE)
        ymin = min(tmp_data$dFF,na.rm=TRUE)


    if(first(unique(tmp_data$protocol)) == "spont" | first(unique(tmp_data$protocol)) == "1Hz") {
            trace_width = 48
            trace_height=12
            xmax_mod = xmax

            
            } else {
                trace_width = 12
                trace_height=8
                xmax_mod = xmax*0.75
                
            }

    print(paste0("ymax is returning the following value :", ymax))

    if (!is.finite(ymax) ) {
        print("Moving onto the next loop I think.")
        next
    }
    if(ymax < 0.2) {
        print("Moving onto next loop because peak is too small for display")
        next
    }
    

    tmp_ID = randomROIs[i]
    .keyvars = rlang::syms(keys)
                    
    tmp_data = tmp_data %>% 
                   mutate(stimKey = paste(!!!.keyvars, sep="-") )

    newPeakColors = list(Signal = "isPeak", Baseline = "noPeak")

    levels(tmp_data$peakColor) <- newPeakColors
    print(levels(tmp_data$peakColor))

    #if()


    tmp_stim_ID = unique(tmp_data$stimKey)  
    stimParadigm_times = def_stim_vlines(tmp_stim_ID)

    vline_switch = !is.null(stimParadigm_times)

    if(add_max == TRUE) {
            usual_frame_time<-median(tmp_data$interFrame,na.rm=TRUE)
            dur_TTL = 0.01 #s

            
            TTL_index = unique(tmp_data$TTL_start)
            max_pt<-tmp_data %>% slice(which.max(dFF)) %>% dplyr::filter(peakID != "NotPeak")
            pct90_line <- as.numeric(0.9*max_pt$dFF)
            pct10_line <- as.numeric(0.1*max_pt$dFF)

            if(is.na(ymax)) {
            
                scalebar_switch = FALSE
                big_scalebar_switch = FALSE

            } else if(ymax > 3 ){

                scalebar_switch = TRUE
                big_scalebar_switch = TRUE

            } else if(ymax > 1.2 ){

                scalebar_switch = TRUE
                big_scalebar_switch = FALSE
             
             } else if(ymax > 0.4) {
                scalebar_switch = FALSE
                big_scalebar_switch = FALSE
                        #} else  if(ymax < 0.4) {
                            #print(tmp_data)
                         #   print(tmp_time_dots)
                          #  next
                        } else { 
                            next
                            }
            
            #desired scalebar size = 1 dFF
            if(big_scalebar_switch == TRUE & scalebar_switch == TRUE){  
                                            scalebar_y_length = 2  
                                            y_length_numeric = scalebar_y_length
                                            #y_length_str = y_length_numeric
                                            add_scalebars = TRUE

                        }    else if(scalebar_switch == TRUE){  
                                            scalebar_y_length = 1  
                                            y_length_numeric = scalebar_y_length
                                            #y_length_str = y_length_numeric
                                            add_scalebars = TRUE
                                             
                                        } else {  
                                            scalebar_y_length = 0.5*ymax
                                            y_length_numeric = round(scalebar_y_length,1)
                                            #y_length_str = y_length_numeric#as.character(y_length_numeric)
                                            add_scalebars = FALSE

                                              }
            scalebar_x_length = 0.1
            x_length_numeric = scalebar_x_length*1000
#            x_length_str = x_length_numeric#as.character(x_length_numeric)


            offset_x = scalebar_x_length*0.052
            user_x = 1.065 #(xmax_mod-1) * 0.1 + 1
            user_y = (ymax*1.2) *0.35

            scalebar_x = data.frame(x_start = user_x,
                                        x_end = user_x+scalebar_x_length,
                                        y_start = user_y,
                                        y_end = user_y)
             scalebar_y =  data.frame(x_start = user_x+offset_x,
                                        x_end = user_x+offset_x,
                                        y_start = user_y,
                                        y_end = user_y+scalebar_y_length)

             label_x = data.frame(#label = bquote(x_length_numeric ~ "ms"),
                                    x = user_x+(scalebar_x_length/2),
                                    y=  user_y)

             label_y = data.frame(x = user_x,
                                    y=  user_y+(scalebar_y_length)/2)

             label_x_expr = list(bquote(.(x_length_numeric)~"ms") )
             label_y_expr = list(bquote(.(y_length_numeric)~Delta*"F/F") )

    }


                             

        if("decay10_time" %in% names(tmp_time_dots)){
            #print(decay_)
            decay10_lineend = tmp_time_dots$decay10_time
        } else if ("absoluteTime" %in% names(max_pt)) {
            decay10_lineend = max_pt$absoluteTime+scalebar_x_length
        } else {
            add_scalebars = FALSE
        }

        vars_1 = c('.fitted','dFF')
        vars_2 = c('rise10_time','rise10_dFF','rise90_time','rise90_dFF','half_rise_time','half_rise_dFF','half_decay_time','half_decay_dFF','decay90_time','decay90_dFF','decay10_time','decay90_dFF')
        if(!all(vars_1 %in% names(tmp_data)) ) {
            print("Didn't find everything in the main dataframe? Skipping to next trace.")
            next
        } else if(!all(vars_2 %in% names(tmp_time_dots))) {
            print("didn't find everything in the time_dots dataframe? Skipping to next trace.")
            next
        
        } else {
        parse.labels <- function(x) parse(text = x)

        text_size = 9.09


        tmp_tracePlot<-ggplot(tmp_data, aes(x=absoluteTime, y = dFF, colour = peakColor, group =1)) +
            
            
            {if(add_scalebars)geom_segment(data=scalebar_x, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 3, colour="black")}+
            {if(add_scalebars)geom_segment(data=scalebar_y, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 3, colour="black")}+
            {if(add_scalebars)geom_label(aes(x=label_x$x,y=label_x$y, label=label_x_expr), nudge_y = (user_y*-0.15), size= text_size, colour="black",parse=TRUE,label.size=NA,family="sans")}+ #x_label_scalebar
            {if(add_scalebars)geom_label(aes(x=label_y$x,y=label_y$y, label=label_y_expr), nudge_x = (user_x*-0.1), size= text_size, colour="black",parse=TRUE,label.size=NA,family="sans")}+ #y_label_scalebar
            
           
            {if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=half_rise_time,y=half_rise_dFF,xend=half_decay_time,yend=half_decay_dFF, colour="t[1/2]"), size= 2, alpha=1)}+ #segment to join the two
            geom_path(size = 1.5,lineend='round',alpha=1.0) + 
            geom_path(aes(x=absoluteTime, y = .fitted, group=1, colour="Exp.~Decay~Fit"),linetype="solid", size=2,alpha=0.8)+
            {if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "grey21", alpha = 1, linetype = "longdash", size =0.9)} +
            {if(vline_switch)geom_label(aes(x = stimParadigm_times, y = ifelse(is.finite(ymax), ymax*1.24, 0 ) ),label = "Stimulus", colour = "black", family='sans', size=text_size,label.size=NA)}+
            
             #half_width
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = half_rise_time, y = half_rise_dFF, colour="t[1/2]"), shape=21, fill="aquamarine4",size= 5, alpha=0.85)}+ #half_t[rise]
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = half_decay_time, y = half_decay_dFF, colour="t[1/2]"), shape=21, fill="aquamarine4",size= 5, alpha=0.85)}+ #half_decay_time


             ##dashed lines indicating 90 and 10% dFF levels
           
            {if(add_scalebars)geom_segment(aes(x=max_pt$absoluteTime-(0.5*scalebar_x_length), xend=max_pt$absoluteTime+(scalebar_x_length), y=pct90_line,yend=pct90_line), colour = "grey41", alpha = 1, linetype = "dashed", size = 0.6)}+
            {if(add_scalebars)geom_segment(aes(x=max_pt$absoluteTime-(0.5*scalebar_x_length), xend=decay10_lineend+(scalebar_x_length), y=pct10_line,yend=pct10_line), colour = "grey41", alpha = 1, linetype = "dashed", size = 0.6)}+

            {if(add_scalebars)geom_label(aes(x = max_pt$absoluteTime+(scalebar_x_length), y = pct90_line), nudge_x = 0.05, label = "90%",fontface="italic", colour = "black", family='sans', size=text_size,label.size=NA,alpha=0.95,fill=NA)}+
            {if(add_scalebars)geom_label(aes(x = tmp_time_dots$decay10_time+(scalebar_x_length), y = pct10_line), nudge_x = 0.05,label = "10%",fontface="italic", colour = "black", family='sans', size=text_size,label.size=NA,alpha=0.95,fill=NA)}+

            
            # #rise_time
            # {if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=rise10_time,y=rise10_dFF,xend=rise90_time,yend=rise90_dFF), size= 3, alpha=0.75,colour="green")}+ #segment to join the two
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = rise10_time, y = rise10_dFF, colour="t[rise]"), shape=21, fill="deepskyblue",size= 5, alpha=0.85)}+ #10riseTime
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = rise90_time, y = rise90_dFF, colour="t[rise]"), shape=21, fill="deepskyblue",size= 5, alpha=0.85)}+ #90riseTime

            #  {if(add_time_dots)geom_segment(data=tmp_time_dots, aes(x=decay10_time,y=decay10_dFF,xend=decay90_time,yend=decay90_dFF), size= 3, alpha=0.75,colour="orange")}+ #segment to join the two
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = decay10_time, y = decay10_dFF, colour="t[decay]"), shape=21, fill="darkorchid", size= 5, alpha=0.85)}+ #10decayTime
            {if(add_time_dots)geom_point(data=tmp_time_dots, aes(x = decay90_time, y = decay90_dFF,colour="t[decay]"), shape=21, fill="darkorchid",size= 5, alpha=0.85)}+ #90decayTime


            #decay_time

            
           
         
            {if(add_max)geom_point(data=max_pt, aes(x=absoluteTime,y=dFF, colour = "Peak[Delta*F/F]"), shape=23, fill="blue",alpha=0.65, size=8)}+
                                  
           
           
            labs( x="Time (s)",
                  y=expression(Delta*"F/F"),
                  title="",#unique(tmp_data$ROI_key),
                  subtitle= "")+#ifelse(!is.null(stimParadigm_times), paste0("First stim estimated at:", round(stimParadigm_times[1],3), " seconds"), NA))+
            scale_colour_manual(breaks = c("Signal", "Baseline","Exp.~Decay~Fit", "t[rise]", "t[decay]", "t[1/2]", "Peak[Delta*F/F]"), 
                                values = c( "Signal" ="black","Baseline" = "grey72","Exp.~Decay~Fit"="red", "t[rise]" = "deepskyblue", "t[decay]" = "darkorchid","t[1/2]" = "aquamarine4", "Peak[Delta*F/F]" = "blue"),
                                labels = parse.labels)+
            #coord_cartesian(xlim=c(0,xmax))+
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
            coord_cartesian(ylim=c(ymin,ymax*1.2),xlim=c(0.75,xmax_mod))+
            #scale_x_continuous(limits=c(1.0,2.0))+
            my.theme+
            {if(title_switch)theme(plot.title=element_text(size=28, family="sans"),
                                    legend.position = c(0.8,0.8),
                                    legend.text.align = 0,
                                        plot.subtitle=element_blank(),
                                        strip.text=element_blank(),
                                        plot.margin = margin(1,1,1,1, "cm"),
                                        axis.line = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.title = element_blank(),
                                        axis.ticks = element_blank())
                                        }+
           
            facet_wrap(~ROI_key, scales ="free",ncol=1,nrow=1)
      
    #save_plot(filename = paste0("traceFits_",i), plot=tmp_tracePlot, plot_width=trace_width*1.5,plot_height=trace_height*1.5, scale_f = scalefactor,dpi_val = 300)
            
   ggsave(filename=paste0("traceFits_",i,".jpeg"),plot=tmp_tracePlot, device="jpeg",bg='white',dpi=dpi_val, units="in",width=trace_width*1.1,height=trace_height*1.1)   
    }
  
   



    rm(list=ls(pattern="tmp_"))
}

toc()

rm(i,n,ROIs,randomROIs,ROIs_to_sample,plotChars,xmax,ymin,ymax,vline_switch,title_switch,trace_width,trace_height)



#warnings back on
options(warn = oldw)

