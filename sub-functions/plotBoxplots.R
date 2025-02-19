subDir <- "peak_violinPlots"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 


library(ggforce)




plotBoxplots = function(df,groupers,plotBy,var, secondAxis=NULL, color_override=NULL) {

			
      color_switch = !is.null(color_override)
			getVar<-df[ ( names(df) %in% var ) ]
      
      #getVar.finite = as.numeric( getVar[[1]] )
      #getVar.finite = getVar.finite[is.finite(getVar.finite)]
      #fix infinite

			#median = median(as.numeric(getVar[[1]]), na.rm=TRUE)
			
      if(var == "releaseProbability") {
        lower_bound_y = 0
        jitter_switch = FALSE
        smooth_switch = TRUE

      } else {
        lower_bound_y = as.numeric( min( getVar ) )
        jitter_switch = FALSE
      }
      upper_bound_y = as.numeric( max( getVar ) )
			dot_y_range = upper_bound_y/50
      #binwidth = upper_bound_y / 30  
			getTitle = rlang::parse_expr(plotBy)	
			getPlotTitle = rlang::parse_expr(var)
			count = "..count.."
			ncount = "..ncount.."

			units = case_when(var == "halfwidth_ms" ~ expression("FWHM (ms)"),
								var == "amplitude" ~ expression(Delta*"F/F"),
								var == "tau_decay_ms" ~ expression(tau[decay]*" (ms)"),
	             	var == "interSpike_ms" ~ expression(Delta*"t (ms)"),
                var == 'releaseProbability' ~ expression(p[syn]),
                var == 'cum_amplitudes' ~ expression('Cumulative '*Delta*'F/F'),
                var == 'n_sync' ~ expression(n[sync]),
                var == 'n_async' ~ expression(n[async]),
                var == 'n_tot' ~ expression(n[total]),
                var == 'async_vs_sync' ~ expression(n[async]/n[sync]),
                var == 'CV' ~ expression("Coefficient of Variance per ROI"),					
                var == 'transmission_rate' ~expression("Synaptic Transmission Rate"),
                var == 'amplitudeGlu_failure' ~expression("GluSnFR3 "*Delta*"F/F"[failure]),
                var == 'amplitudeGlu_success' ~expression("GluSnFR3 "*Delta*"F/F"[success]),
                var == 'amplitudeJF_success'~expression("JF646 "*Delta*"F/F"[success]),
                var == 'interSpikeCheck' ~expression(Delta*t[success]) ) 
            
			

      if(missing(secondAxis)) {

			boxPlot <-ggplot(df, aes_string(x=plotBy,y=var, colour=plotBy)) +
                            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75),fill="white",size=1.5, alpha=0.4)+ 
                            geom_sina(alpha=0.2)+
                            # {if(jitter_switch==FALSE)geom_dotplot(aes_string(colour=plotBy, fill=plotBy),
                            #                                       binaxis= "y",
                            #                                        stackdir = "center",
                            #                                        dotsize = 0.5,
                            #                                        alpha=0.5,
                            #                                        binwidth=dot_y_range,
                            #                                        colour=NA) } +
                            
                            labs( x="",
                                  y=units,
                                  title=bquote(.(getPlotTitle)),
                                  subtitle="")+
                            coord_cartesian(ylim = c(lower_bound_y,upper_bound_y))+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(values=color_override)}+
                            {if(color_switch)scale_fill_manual(values =color_override)}+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")#+
                            #facet_grid( eval( expr( ~!!ensym(secondAxis) ) ) )

               			ggsave(filename=paste(getPlotTitle,getTitle,"_boxplot.png", sep='_'),plot=boxPlot, device="png",dpi=600, units="in",width=14,height=16)


        } else { 
          #print("boxplot has secondAxis")
          getSecondAxis = rlang::parse_expr(secondAxis)
          boxPlot <-ggplot(df, aes_string(x=plotBy,y=var, colour=plotBy)) +
                            geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill="white",size=1.5, alpha=0.4)+ 
                            geom_sina(alpha=0.2)+
                            # {if(jitter_switch==FALSE)geom_dotplot(aes_string(fill=plotBy),
                            #                                       binaxis= "y",
                            #                                        stackdir = "center",
                            #                                        dotsize = 0.5,
                            #                                        alpha=0.5,
                            #                                        binwidth=dot_y_range,
                            #                                        colour=NA) } +
                            
                            labs( x="",
                                  y=units,
                                  title=bquote(.(getPlotTitle)~"by"~.(getSecondAxis)),
                                  subtitle="")+
                            coord_cartesian(ylim = c(lower_bound_y,upper_bound_y))+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(values=color_override)}+
                            {if(color_switch)scale_fill_manual(values =color_override)}+
                            theme_tufte()+
                            my.theme+
                            theme(legend.position = "none")+
                            facet_grid( eval( expr( ~!!ensym(secondAxis) ) ) )

                    ggsave(filename=paste(getPlotTitle,getTitle,getSecondAxis,"_violinplot.png", sep='_'),plot=boxPlot, device="png",dpi=600, units="in",width=14,height=16)



        }

}
