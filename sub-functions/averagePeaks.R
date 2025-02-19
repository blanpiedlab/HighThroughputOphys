#averagePeaks.R

#function(df,groupers,plotBy){}
#accepts a dataframe and outputs faceted, averaged peaks


subDir <- "averagedPeaks"

source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 




averagePeaks<- function(df,groupers,plotBy,levels, secondAxis=NULL, color_override = NULL, thirdAxis = NULL){

            interFrame = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)

            color_switch = !is.null(color_override)
            facet_switch = !is.null(secondAxis) 
      
            justPeaks<-df %>%
                        group_by_at(groupers) %>%
                        dplyr::filter(peakID != "NotPeak", fitQualityPass == TRUE) %>%
                        mutate(normTime = absoluteTime - min(absoluteTime))



            justPeaks[[plotBy]]<- factor(justPeaks[[plotBy]],
                               levels = levels)

            suppressMessages(getMedian_peakLength<- justPeaks %>% group_by_at(groupers) %>% summarise(max_duration = max(normTime)))


            upper_bound_x <- max(justPeaks$normTime)
            
            frameRate<-seq(0,upper_bound_x, by=interFrame)
            numBins = ceiling( length(frameRate) * 1.5)
            
            upper_bound_x_exclusive<-median(getMedian_peakLength$max_duration)
            frameRate_exclusive<-seq(0,upper_bound_x_exclusive, by=interFrame)
            numBins_exclusive = ceiling( length(frameRate_exclusive) * 1.5) 
            getTitle = rlang::parse_expr(plotBy)


            averagePeakPlot<-ggplot(justPeaks, aes(x=normTime, y=dFF))+
                            geom_line(aes(group=ROIpeaks), alpha=0.7,colour='grey',size=1)+
                            {if(facet_switch == FALSE)stat_summary_bin(aes_string(group = plotBy, colour=plotBy), geom="smooth", fun = mean, bins=numBins_exclusive)}+
                            {if(facet_switch)stat_summary_bin(aes_string(group = interaction(plotBy,secondAxis), colour=plotBy), geom="smooth", fun = mean, bins=numBins_exclusive)}+
                            labs( x="Normalized Time (s)",
                                    y=expression(Delta*"F/F"),
                                    title=bquote("Comparison of average peaks by :"~.(getTitle))
                                    )+
                            
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(values=color_override)}+
                            {if(facet_switch == FALSE)facet_grid(eval(expr(~!!ensym(plotBy))))}+
                            {if(facet_switch)facet_grid(eval(expr(!!ensym(secondAxis)~!!ensym(plotBy))))}+
                            theme_tufte()+
                            my.theme+
                            coord_cartesian(xlim=c(0,upper_bound_x))+
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

            ggsave(filename=paste0("_averagedPeaks.png"),plot=averagePeakPlot, device="png",dpi=600, units="in",width=20,height=12)


            averagePeakPlot<-ggplot(justPeaks, aes(x=normTime, y=dFF))+
                            geom_line(aes(group=ROIpeaks), alpha=0.7,colour='grey',size=1)+
                            {if(facet_switch == FALSE)stat_summary_bin(aes_string(group = plotBy, colour=plotBy), geom="smooth", fun = mean, bins=numBins_exclusive)}+
                            {if(facet_switch)stat_summary_bin(aes_string(group = interaction(plotBy,secondAxis), colour=plotBy), geom="smooth", fun = mean, bins=numBins_exclusive)}+
                            labs( x="Normalized Time (s)",
                                    y=expression(Delta*"F/F"),
                                    title=bquote("Excluding peaks greater than the median peak duration, comparison of average peaks by :"~.(getTitle))
                                    )+
                            
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(values=color_override)}+
                            {if(facet_switch == FALSE)facet_grid(eval(expr(~!!ensym(plotBy))))}+
                            {if(facet_switch)facet_grid(eval(expr(!!ensym(secondAxis)~!!ensym(plotBy))))}+
                            theme_tufte()+
                            my.theme+
                            coord_cartesian(xlim=c(0,upper_bound_x_exclusive))+
                            scale_x_continuous(limits = c(0,upper_bound_x_exclusive))+
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

            ggsave(filename=paste0("exclusive_averagedPeaks.png"),plot=averagePeakPlot, device="png",dpi=600, units="in",width=20,height=12)






}



