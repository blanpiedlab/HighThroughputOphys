#averagePeaks_PPedit.R

#function(df,groupers,plotBy){}
#accepts a dataframe and outputs faceted, averaged peaks

#new edit should include a mechanism by which to facet by an interaction
#Want to look at entire traces via an adaptive grouping... inclusive of a cutoff that makes sense




subDir <- "averagedPeaks_PPedit"

source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 




averagePeaks<- function(df,groupers,plotBy,levels, secondAxis=NULL, color_override = NULL,  thirdAxis = NULL,thirdAxis_levels = NULL){

            interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)

            cumplus <- function(y) Reduce(function(a,b) a + b > 0, y, 0, accum=TRUE)[-1]

            
            color_switch = !is.null(color_override)
            facet_switch = !is.null(secondAxis)
            thirdAxis_switch = !is.null(thirdAxis) 
            #interaction_swtich = !is.null(interactors)        
      
            justPPTrace<-df %>%
                        group_by_at(groupers) %>%
                        #mutate( PPstart = case_when(index == 180  ~ 1,
                         #                               TRUE  ~ 0),
                          #      PPexit = 0,
#
 #                               cumplus = cumplus(PPstart - PPexit),
  #                              temp = cumplus - c(0,pmin(0,diff(cumplus)))
   #                 
    #                            ) %>%
                        dplyr::filter(index >= 180)

            drops <- c("PPstart","PPexit","temp")  #extraneous data.frame columns clogging up memory
            justPPTrace<- justPPTrace[ , !(names(justPPTrace) %in% drops)]

            justPPTrace<- justPPTrace %>% mutate(neuronCa = paste(neuronSegment,Ca, sep='-'))
            justPPTrace$neuronCa<- factor(justPPTrace$neuronCa, levels = thirdAxis_levels)

                                    


            justPPTrace[[plotBy]]<- factor(justPPTrace[[plotBy]],
                               levels = levels)

            #suppressMessages(getMedian_peakLength<- justPeaks %>% group_by_at(groupers) %>% summarise(max_duration = max(normTime)))


            upper_bound_x <- max(justPPTrace$absoluteTime)
            lower_bound_x <- min(justPPTrace$absoluteTime[which(justPPTrace$index == 180)])
            upper_bound_y <-0.60*max(justPPTrace$dFF)
            lower_bound_y <-min(justPPTrace$dFF)
            frameRate<-seq(lower_bound_x,upper_bound_x, by=interFrame.var) 
           numBins = ceiling( length(frameRate) * 0.5)
            
            #upper_bound_x_exclusive<-median(getMedian_peakLength$max_duration)
            #frameRate_exclusive<-seq(0,upper_bound_x_exclusive, by=interFrame)
            #numBins_exclusive = ceiling( length(frameRate_exclusive) * 1.5) 
            getTitle = rlang::parse_expr(plotBy)

            baseHeight = 4
            baseWidth = 6
            scaledHeight = baseHeight*0.60*length(unique(justPPTrace$neuronCa) )
            scaledWidth = baseWidth*0.60*length( unique(justPPTrace$protocol) )
            text_scale = ceiling(scaledWidth*scaledHeight/20)

            averagePeakPlot<-ggplot(justPPTrace, aes(x=absoluteTime, y=dFF))+
                            geom_line(aes(group=ROIpeaks), alpha=0.7,colour='grey',size=1)+
                            {if(thirdAxis_switch)stat_summary_bin(aes(group=neuronCa, colour=neuronCa), geom="smooth", fun = mean, bins=numBins)}+
                            #{if(thirdAxis_switch)stat_summary_bin(aes_string(group = interaction(plotBy,secondAxis), colour=plotBy), geom="smooth", fun = mean, bins=numBins_exclusive)}+
                            labs( x="Normalized Time (s)",
                                    y=expression(Delta*"F/F"),
                                    title=bquote("Comparison of average peaks by :"~.(getTitle))
                                    )+
                            
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            #{if(color_switch)scale_colour_manual(values=color_override)}+
                            #{if(facet_switch == FALSE)facet_grid(eval(expr(~!!ensym(plotBy))))}+
                            {if(thirdAxis_switch)facet_grid(neuronCa~protocol) }+
                            theme_tufte()+
                            my.theme+
                            coord_cartesian(xlim=c(lower_bound_x,upper_bound_x),ylim=c(lower_bound_y,upper_bound_y))+
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
                                          strip.text = element_text(colour="black", size = 16, family = "Serif")
                                          )

            ggsave(filename=paste0("_averagedPeaks.png"),plot=averagePeakPlot, device="png",dpi=600, units="in",width=scaledWidth,height=scaledHeight)


          

}



