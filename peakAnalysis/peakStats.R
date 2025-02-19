#peakStats.R




peakStats<- function(df,groupers,subgroup,plotBy, varList, levels, secondAxis = NULL, color_override=color_override){
						
						secondAxis_switch = !is.null(secondAxis) #check if secondAxis IS NOT NULL. If secondAxis exists, this secondAxis_switch will return TRUE
						.group_vars = ifelse(secondAxis_switch == TRUE, c(plotBy, secondAxis), plotBy) #we'll come back to this

						
						 taus<-suppressMessages(df %>% group_by_at(groupers) %>%
					                 unnest(tidied) %>%
					                 spread(term,estimate) %>%
					                 ungroup() %>%
					                 group_by_at(groupers) %>%
					                 summarise(	tau_decay = unique( tau_decay[!is.na(tau_decay)] )																				#tau_decays that are not NA, one for each peak. 
					                            )
					                 )

					   	onlyCleanPeaks<- suppressMessages(df %>% group_by_at(groupers) %>% dplyr::filter(#fitQualityPass == TRUE, 
					   																					peakID != "NotPeak") %>%#, stimEpoch == 1) %>%
					   																		#dplyr::filter(absoluteTime <= ( min(absoluteTime)+0.150 ) )   %>% ##from here, we are only looking at times after the first stimulus. only grab times <500 ms since stimulus.
                                                                            
					   										summarise(amplitude = max(dFF),
					   													interSpike = unique(interSpike)#,
					   													#maxTime = absoluteTime[which(dFF==amplitude)],
					   													#firstTime = maxTime-interSpike 
					   										)
					   									) 
					   	onlyCleanPeaks$stimEpoch[which(is.na(onlyCleanPeaks$stimEpoch))]<- 0 #coerce NAs might be necessary? 

#					   	onlyCleanPeaks
						
						# colnames_toCheck = C("maxTime","riseTime") 

						# if(colnames_toCheck %in% colnames(onlyCleanPeaks)) {
						# 	print("maxTime and riseTime seem to be in the dataset.")
						# } 
						
					   	peakStats<- suppressMessages(left_join(onlyCleanPeaks, taus) %>% group_by_at(groupers) %>% 
					   											summarise(tau_decay_ms = tau_decay*1000,
					   														interSpike_ms = interSpike*1000,
					   														amplitude = amplitude#,
					   														#riseTime_ms = riseTime*1000
					   														) %>%
					   											ungroup() %>%
					   											group_by_at(subgroup) %>% 
					   											 dplyr::filter(	tau_decay_ms > 10, 
					   											  				tau_decay_ms < quantile(tau_decay_ms, 0.95, na.rm=TRUE), 
																			  	tau_decay_ms > quantile(tau_decay_ms, 0.05, na.rm=TRUE),
																			  	tau_decay_ms < 150,

																			  	amplitude < quantile(amplitude, 0.95, na.rm=TRUE), 
																			  	amplitude > quantile(amplitude, 0.05, na.rm=TRUE),

																			  	interSpike_ms < 200) 
					   									)



						
						peakStats[[plotBy]]<- factor(peakStats[[plotBy]],
                               									levels = levels)

						
						 

						 # averagePeakPlot<-ggplot(traces, aes(x=normTime, y=dFF))+
       #                      #geom_ribbon(aes(x = absoluteTime, ymax = Inf, fill = category), ymin=0, alpha=0.3,size=0) +                #uncomment this if you want to visualize the stim Epochs as a color-coded ribbon
       #                      #geom_line(aes(group=ROI_keys), alpha=0.4,colour='grey',size=1.5)+
       #                      {if(facet_switch)stat_summary_bin(aes_string(group=plotBy, colour=plotBy), geom="smooth", fun = mean, bins=numBins,size=3)}+
       #                      geom_vline(xintercept = 0, lty='dashed',colour='black',size=2)+
                            
       #                      labs( x="Normalized Time (s)",
       #                              y=expression(Delta*"F/F"),
       #                              title="" #bquote("Comparison of average peaks by :"~.(getTitle))
       #                              )+
                            
       #                      {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
       #                      #scale_fill_manual(values = getPalette(colourCount))+
       #                      {if(color_switch)scale_colour_manual(labels=Ca_labels,values=color_override)}+
       #                      #{if(facet_switch == FALSE)facet_grid(eval(expr(~!!ensym(plotBy))))}+
       #                      #facet_grid(eval( expr( !!ensym(secondAxis) ~ !!ensym(plotBy) ) ), scales = 'free_x')+
       #                      theme_tufte()+
       #                      my.theme+
       #                      scale_x_continuous(labels = number_format(accuracy = 0.01), breaks=c(0.00,0.25,0.50))+
       #                      scale_y_continuous(labels = number_format(accuracy = 0.1), breaks =c(0.0,0.5,1.0,1.5,2.0))+
       #                      theme(legend.position = c(.95, .95),
       #                              legend.justification = c("right", "top"),
       #                              legend.box.just = "right",
       #                              legend.margin = margin(6, 6, 6, 6)
       #                              )+
       #                      coord_cartesian(ylim=c(lower_bound_y,2.1), xlim=c(lower_bound_x,0.6))
                           
       #   			   ggsave(filename=paste0("_averagedPeaks_nofacets_.png"),plot=averagePeakPlot, device="png",dpi=600, units="in",width=16,height=12)



						

						if(secondAxis_switch == FALSE) {
							

							print("No second axis specified.")
							

							source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotBoxplots.R")

							
							for (i in 1:length(varList)) { 

							var = varList[i]
							plotBoxplots(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, color_override=color_override)

							}	




							source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotHistos.R") 
							
							
							for (j in 1:length(varList)) { 

							var = varList[j]
							plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var,color_override=color_override)

							}	


						} else {

							print("Plotting graphs by a second variable.")
							print(secondAxis)

							
							source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotBoxplots.R")

							
							for (i in 1:length(varList)) { 

							var = varList[i]
							plotBoxplots(df = peakStats, groupers = subgroup ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

							}	




							source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotHistos.R") 
							
							
							for (j in 1:length(varList)) { 

							var = varList[j]
							plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

							}	

							

						

						}


							
						

						
						peakStats


}











