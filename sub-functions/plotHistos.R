


subDir <- "peakHistos"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 


library(ggforce)
library(scales)







plotHistos = function(df,groupers,plotBy,var, secondAxis = NULL, color_override = NULL, bin_override = NULL,title_string=NULL) {

			Ca_labels = c(expression("0.5 mM "*Ca^'2+'),#expression("1 mM "*Ca^'2+'), 
        expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+'), expression("8 mM "*Ca^'2+') )


      plot_height = 16
      plot_width = 16

      color_switch = !is.null(color_override)
      bin_switch = !is.null(bin_override) 
      
      title_switch = !is.null(title_string) 
      

      if(missing(secondAxis)) {
          measureVars<- c(plotBy,var)
        }else{
          secondAxis=secondAxis
          measureVars<- c(plotBy,secondAxis,var) 
        }
			getVar<-df[ ( names(df) %in% measureVars ) ]
      varRange = df[ ( names(df) %in% var ) ] 
      lower_bound_x = as.numeric( min( as.numeric( varRange[[1]] ) ) )
      upper_bound_x = as.numeric( max( as.numeric( varRange[[1]] ) ) )
      



			if(missing(secondAxis)) {
          grouping_vars = plotBy
          my_medians = suppressMessages( getVar %>% group_by_at(.vars=grouping_vars) %>% na.omit() %>% summarise(median_val = median(!!ensym(var),na.rm=TRUE) ) %>% mutate(labels = round(median_val,2),x.pos = 0.5*(upper_bound_x-median_val)+median_val ) ) 
          
      
        }else{
          #print("Attempting to set the medians including secondAxis:")
          #print(secondAxis)
          #print("These are the grouping_vars fed to group_by_at()")
          grouping_vars = c(plotBy,secondAxis)
          #print(grouping_vars)
         
          my_medians = suppressMessages( getVar %>% dplyr::group_by_at(.vars=grouping_vars) %>% summarise(median_val = median(!!ensym(var),na.rm=TRUE) ) %>% mutate(labels = round(median_val,2), x.pos = 0.5*(upper_bound_x-median_val)+median_val  ) ) 
         
        }
      



      #medians = getVar %>% group_by_at(plotBy,secondAxis)

      #median = median(as.numeric(getVar[[1]]), na.rm=TRUE)
			binwidth = (upper_bound_x - lower_bound_x) / 50  
			getTitle = rlang::parse_expr(plotBy)
      getPlotTitle = rlang::parse_expr(var)
			count = "..count.."
			ncount = "..ncount.."

			units = case_when(var == "halfwidth_ms" ~ expression("FWHM (ms)"),
								var == "amplitude" ~ expression(Delta*"F/F"),
								var == "tau_decay_ms" ~ expression(tau[decay]*" (ms)"),
								var == "interSpike_ms" ~ expression(Delta*"t (ms)"),
                var == 'transmission_rate' ~expression("Synaptic Transmission Rate"),
                var == 'amplitudeGlu_failure' ~expression("GluSnFR3 "*Delta*"F/F"[failure]),
                var == 'amplitudeGlu_success' ~expression("GluSnFR3 "*Delta*"F/F"[success]),
                var == 'amplitudeJF_success'~expression("JF646 "*Delta*"F/F"[success]),
                var == 'interSpikeCheck' ~expression(Delta*t[success])  
                )

      
        if(missing(secondAxis)) {
              

                  

                  histoCount <-ggplot(df, aes_string(x=var, colour=plotBy, fill=plotBy)) +
                                        {if(bin_switch == FALSE)geom_histogram(aes(y=..count..),binwidth = binwidth,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                        {if(bin_switch == FALSE)geom_freqpoly(data=df, aes_string(x=var, y=count),binwidth = binwidth, alpha=0.8,,boundary =0,size=2)}+
                                        {if(bin_switch)geom_histogram(aes(y=..count..),binwidth = bin_override,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                        {if(bin_switch)geom_freqpoly(data=df, aes_string(x=var, y=count),binwidth = bin_override, alpha=0.8,,boundary =0,size=2)}+
                                        
                                        geom_vline(data = my_medians, aes(xintercept = median_val), colour="black",lty="dashed",size=2.5)+
                                        # annotate(geom= "text",
                                        #           label=paste0(expression("Median "*Delta*"F/F: "), 
                                        #                    as.character(round(my_medians$median_val,1)),
                                        #                    col='black', size=4, x=(100))+
                                                                            
                                        labs( x=units,
                                              y=expression("Count"),
                                              title="",#bquote("Comparing by :"~.(getTitle)),
                                              subtitle="")+
                                        scale_x_continuous(limits = c(lower_bound_x,upper_bound_x))+
                                        scale_y_continuous(expand = c(0,0),labels = number_format(accuracy = 0.1))+
                                        {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                        {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                        {if(color_switch)scale_colour_manual(values=color_override)}+
                                        {if(color_switch)scale_fill_manual(values=color_override, labels=Ca_labels)}+
                            
                                        theme_tufte()+
                                        my.theme+
                                        theme(legend.position = "none"#,
                                              #panel.spacing = unit(2, "lines")
                                              )+
                                        facet_grid(eval(expr(!!ensym(plotBy)~.)) 
                                )

                                if(title_switch){ ggsave(filename=paste(title_string,getPlotTitle,getTitle,"_countHisto.png", sep='_'),plot=histoCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height) 
                                                } else {
                                ggsave(filename=paste(getPlotTitle,getTitle,"_countHisto.png", sep='_'),plot=histoCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height)
                                                }



                            
                  histoNCount <-ggplot(df, aes_string(x=var, colour=plotBy, fill=plotBy)) +
                                        {if(bin_switch == FALSE)geom_histogram(aes(y=..ncount..),binwidth = binwidth,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                        {if(bin_switch == FALSE)geom_freqpoly(data=df, aes_string(x=var, y=ncount),binwidth = binwidth, alpha=0.8,,boundary =0,size=2)}+
                                        {if(bin_switch)geom_histogram(aes(y=..ncount..),binwidth = bin_override,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                        {if(bin_switch)geom_freqpoly(data=df, aes_string(x=var, y=ncount),binwidth = bin_override, alpha=0.8,,boundary =0,size=2)}+
                                        #geom_histogram(aes(y=..ncount..),binwidth = binwidth,alpha=0.5, boundary = 0,position="identity",colour="black")+
                                        
                                        #geom_freqpoly(data=df, aes_string(x=var, y=ncount),binwidth = binwidth, alpha=0.8,,boundary =0,size=1)+
                                        geom_vline(data = my_medians, aes(xintercept = median_val), colour="black",lty="dashed",size=2.5)+
                                       # geom_text(data = my_medians, aes(label = labels,x =x.pos), col='black', size=8,y=0.9)+



                                        # annotate(geom= "text",
                                        #            label=paste0(expression("Median "*Delta*"F/F: "), 
                                        #                     as.character(round(my_medians$median_val,1))),
                                        #                     col='black', size=4, x=(0.8*(upper_bound_x-lower_bound_x)),y=0.9)+
                                                                            
                                        labs( x=units,
                                              y=expression("N/N"[max]),
                                              title="",#bquote("Comparing by :"~.(getTitle)),
                                              subtitle="")+
                                        scale_x_continuous(limits = c(lower_bound_x,upper_bound_x),labels = number_format(accuracy = 0.1))+
                                        scale_y_continuous(expand = c(0,0),breaks=c(0,0.5,1.0),labels = number_format(accuracy = 0.1))+
                                        {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                        {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                        {if(color_switch)scale_colour_manual(values=color_override)}+
                                        {if(color_switch)scale_fill_manual(values=color_override, labels=Ca_labels)}+
                                        guides(colour="none")+
                                        theme_tufte()+
                                        my.theme+
                                        theme(strip.text=element_blank(),legend.position = "right"#,
                                              #panel.spacing = unit(2, "lines")
                                              )+
                                        facet_grid(eval(expr(!!ensym(plotBy)~.))  
                                )

                                
                                if(title_switch){ ggsave(filename=paste(title_string,getPlotTitle,getTitle,"_NcountHisto.png", sep='_'),plot=histoNCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height) 
                                                } else {
                                ggsave(filename=paste(getPlotTitle,getTitle,"_NcountHisto.png", sep='_'),plot=histoNCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height)
                                                }




                                
                          

              } else {
                        #print("histo has second axis")
                        getSecondAxis = rlang::parse_expr(secondAxis) 
                        
                        histoCount <-ggplot(df, aes_string(x=var, colour=plotBy, fill=plotBy)) +
                                                  {if(bin_switch == FALSE)geom_histogram(aes(y=..count..),binwidth = binwidth,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                                  {if(bin_switch == FALSE)geom_freqpoly(data=df, aes_string(x=var, y=count),binwidth = binwidth, alpha=0.8,,boundary =0,size=2)}+
                                                  {if(bin_switch)geom_histogram(aes(y=..count..),binwidth = bin_override,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                                  {if(bin_switch)geom_freqpoly(data=df, aes_string(x=var, y=count),binwidth = bin_override, alpha=0.8,,boundary =0,size=2)}+
                                                  
                                                  geom_vline(data=my_medians,aes(xintercept = median_val), colour="black",lty="dashed",size=2.5)+
                                                  
                                                  labs( x=units,
                                                        y=expression("Count"),
                                                        title="",#bquote("Comparing by :"~.(getTitle)~"vs."~.(getSecondAxis)),
                                                        subtitle="")+
                                                  scale_x_continuous(limits = c(lower_bound_x,upper_bound_x))+
                                                  scale_y_continuous(expand = c(0,0),labels = number_format(accuracy = 0.1))+
                                                  {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                                  {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                                  {if(color_switch)scale_colour_manual(values=color_override)}+
                                                  {if(color_switch)scale_fill_manual(values=color_override, labels=Ca_labels)}+
                                                  theme_tufte()+
                                                  my.theme+
                                                  theme(legend.position = "none"#,
                                                        #panel.spacing = unit(2, "lines")
                                                        )+
                                                  facet_grid(eval(expr(!!ensym(plotBy)~!!ensym(secondAxis)))
                                          )



                                             if(title_switch){ ggsave(filename=paste(title_string,getPlotTitle,getTitle,getSecondAxis,"_countHisto.png", sep='_'),plot=histoCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height) 
                                                } else {
                                                              ggsave(filename=paste(getPlotTitle,getTitle,getSecondAxis,"_countHisto.png", sep='_'),plot=histoCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height)
                                                }




                                          




                                      
                            histoNCount <-ggplot(df, aes_string(x=var, colour=plotBy, fill=plotBy)) +
                                                  {if(bin_switch == FALSE)geom_histogram(aes(y=..ncount..),binwidth = binwidth,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                                  {if(bin_switch == FALSE)geom_freqpoly(data=df, aes_string(x=var, y=ncount),binwidth = binwidth, alpha=0.8,,boundary =0,size=2)}+
                                                  {if(bin_switch)geom_histogram(aes(y=..ncount..),binwidth = bin_override,alpha=0.5, boundary = 0,position="identity",colour="black",size=1.5)}+
                                                  {if(bin_switch)geom_freqpoly(data=df, aes_string(x=var, y=ncount),binwidth = bin_override, alpha=0.8,,boundary =0,size=2)}+
                                                  geom_vline(data=my_medians,aes(xintercept = median_val), colour="black",lty="dashed",size=2.5)+
                                                  #geom_text(data = my_medians, aes(label = labels, x=x.pos), col='black', size=8,y=0.9)+

                                                  # annotate(geom= "text",
                                                  #           label=paste0(expression("Median "*Delta*"F/F: "), 
                                                  #                 as.character(round(my_medians$median_val,1))),
                                                  #           col='black', size=4, x=(0.8*(upper_bound_x-lower_bound_x)),y=0.9)+
                                        
                                                  labs( x=units,
                                                        y=expression("N/N"[max]),
                                                        title="",#bquote("Comparing by :"~.(getTitle)~"vs."~.(getSecondAxis)),
                                                        subtitle="")+
                                                  scale_x_continuous(limits = c(lower_bound_x,upper_bound_x),labels = number_format(accuracy = 0.1))+
                                                  scale_y_continuous(expand = c(0,0),breaks=c(0,0.5,1.0),labels = number_format(accuracy = 0.1))+
                                                  {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                                  {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                                  {if(color_switch)scale_colour_manual(values=color_override)}+
                                                  {if(color_switch)scale_fill_manual(values=color_override, labels=Ca_labels)}+
                                                  guides(colour="none")+
                                                  theme_tufte()+
                                                  my.theme+
                                                  theme(strip.text=element_blank()#,#legend.position = "none",
                                                        #panel.spacing = unit(2, "lines")
                                                        )+
                                                  facet_grid(eval(expr(!!ensym(plotBy)~!!ensym(secondAxis))) 
                                          )





                                             if(title_switch){ ggsave(filename=paste(title_string,getPlotTitle,getTitle,getSecondAxis,"_NcountHisto.png", sep='_'),plot=histoNCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height) 
                                                } else {
                                                              ggsave(filename=paste(getPlotTitle,getTitle,getSecondAxis,"_NcountHisto.png", sep='_'),plot=histoNCount, device="png",dpi=600, units="in",width=plot_width,height=plot_height)
                                                }


                                          
                            
              }



histoNCount
    

}
