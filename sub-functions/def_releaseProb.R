#def_releaseProb.R

#def_releaseProb.R is a script that reads in fitsCleaned, the dataframe with all of the peaks that pass inspection. 
#This should be cross-referenced with peakStats, which separates out the 5-95 quantiles (further filtering).

#The desired output is a set of plots and a dataframe encapsulating the interesting, per-stimulus and/or per-ROI information. 

#Step 1. Read in DF. 
#Step 2. identify release probability on a per-stimulus basis. 
#Step 3. start parallel operations
#Step 3a. release probability on per-stimulus basis graph
#step 3b. cumulative amplitude on per-stimulus basis + average

#step 4. Total N
#step 4a. Nasync, Nsync, Ntot on a per ROI basis - look across stimulus paradigms - this will be corrupted without proper stimulus timing
#step 4b. Total release probability on a per ROI basis (summary stat)



def_releaseProb = function(df_all=sync_async_fitsLabeled_stimKeyed,df_clean=fitsCleaned, 
														groupers_byPeak=groupers_byPeak, 
														groupers_byEpoch=groupers_byEpoch,
														groupers_byROI = groupers_byROI,
														keysROI=NULL,
														keysEpoch=NULL,
														plotBy = plotBy,
														levels = levels,
														varList = varList,
														secondAxis = NULL,
														all_ROIs_override = NULL) {


						all_ROIs_switch = !is.null(all_ROIs_override)
						facet_switch = !is.null(secondAxis)

#						.keyvars_ROI = rlang::syms(keysROI)
#						.keyvars_epoch = rlang::syms(keysEpoch)

						
label0= "peakID"
label1= "stimEpoch"

forROI = c(label0,label1)

						#groupers_byPeak = groupersStimuli
						#groupers_byEpoch = setdiff(groupers_byPeak, label0)
						#groupers_byROI = setdiff(groupers_byEpoch, label1)

						amplitude_perStimulus <- suppressMessages(df_clean %>%
																													dplyr::filter(!is.na(stimEpoch)) %>% 
																													group_by_at(groupers_byPeak) %>%
																																	summarise(amplitude = max(dFF, na.rm=TRUE)) 
																																	)




						releaseProb_perStimulus <- suppressMessages(df_clean %>%
																													dplyr::filter(!is.na(stimEpoch)) %>% 
																													group_by_at(groupers_byEpoch) %>%
																																	summarise(releaseProb = n_distinct(peakID),
																																						n_async = n_distinct(peakID[which(timeClass == "async")]),
																																						n_sync = n_distinct(peakID[which(timeClass == "sync")])
																																						) 
																																	)																					
						all_stimulus_epochs <- suppressMessages( df_all %>%
																												dplyr::filter(!is.na(stimEpoch)) %>%
																												group_by_at(groupers_byEpoch) %>%
																												summarise( epoch_keys = unique(epoch_keys) )
																												)

						allEpochs_releaseProb <- suppressMessages( left_join(all_stimulus_epochs, releaseProb_perStimulus) %>%
																						#mutate(ROI_keys = paste(!!!.keyvars_ROI, sep="-")) %>%
																						arrange(ROI_keys,stimEpoch) )

						allEpochs_releaseProb$releaseProb[is.na(allEpochs_releaseProb$releaseProb)] <- 0
																	
						allEpochs<- suppressMessages( left_join(allEpochs_releaseProb,amplitude_perStimulus) )

						allEpochs$amplitude[is.na(allEpochs$amplitude)] <- 0

						allEpochs$amplitude

						allEpochs_final <- allEpochs %>% group_by_at(groupers_byROI) %>% mutate(cum_amplitude = cumsum(amplitude))


						summaryStats <- allEpochs %>% ungroup() %>% group_by_at(groupers_byROI) %>% dplyr::summarise(releaseProbability = mean(releaseProb),
																																																																								mean_amplitudes = mean(amplitude, na.rm=TRUE),
																																																																								sd_amplitudes = sd(amplitude, na.rm=TRUE),
																																																																								CV = sd_amplitudes/mean_amplitudes,
																																																																								cum_amplitudes = sum(amplitude,na.rm=TRUE),
																																																																								n_async = sum(n_async,na.rm=TRUE),
																																																																								n_sync = sum(n_sync,na.rm=TRUE),
																																																																								n_tot = sum(releaseProb,na.rm=TRUE),
																																																																								async_vs_sync = n_async/n_sync )


						
						#print("No second axis specified.")
							

							source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3 - PTPsigma/plotBoxplots.R")

							#groupers = c('ROI_keys')
						for (i in 1:length(varList)) { 

							var = varList[i]
							plotBoxplots(df = summaryStats, groupers = groupers_byROI,plotBy = plotBy, secondAxis=secondAxis, var = var)

							}	







																	#summarise(amplitude = max(dFF, na.rm=TRUE),
																				#releaseProb = n_distinct(peakID)/1
																	#			) 
																	#		)


																				#%>%
																	#ungroup() %>%
																	#group_by_at(groupersROI) %>%
																	#arrange(stimEpoch) %>%
																	#mutate(ROI_keys = paste(!!!.keyvars, sep="-"),
																#					cum_amplitude = cumsum(amplitude)) #%>%

																					#releaseProb = n_distinct(peakID)/1
																						#)

																			

						#releaseStats_perROI <- suppressMessages(releaseStats_perStimulus %>% group_by_at(groupersROI) %>%
						#																					summarise(releaseProb_mean = mean(releaseProb,na.rm=TRUE),
						#																								cum_amplitude = cumsum(amplitude),
						#																								n_tot = row_number(),
						#																								n_async = row_number(timeClass_keep == "async"),
						#																								n_sync = row_number(timeClass_keep == "sync")
						#																								)
						#																							) 


						#Cumulative plots



#boxplots of cumulative dFF per ROI
#boxplots of overall release probability (#events/#stimuli  on a per ROI basis)
#boxplots of n_async /n_tot per ROI


makeGraphs = TRUE


if (makeGraphs == TRUE){
numBins = max(allEpochs_final$stimEpoch)
						cumPlot <-ggplot(allEpochs_final, aes(x=stimEpoch,y=cum_amplitude)) +
														{if(all_ROIs_switch == FALSE)geom_line(aes(group=ROI_keys), alpha=0.5,colour='grey',size=1)}+
														geom_smooth(aes(group=protocol, colour=protocol, fill=protocol), formula = y~x, method="loess", se=TRUE, span=0.25)+
														#stat_summary_bin(aes(group = protocol, colour=protocol), geom="point", fun = mean, bins=numBins)+
														#stat_summary_bin(aes(group = protocol, colour=protocol), geom="line", fun = mean, bins=numBins)+
                        		#stat_summary_bin(aes(group = protocol, colour=protocol), geom="errorbar", fun = mean_se, bins=numBins)+
                        		scale_colour_brewer(palette = "Dark2")+
                            scale_fill_brewer(palette="Dark2")+
                            labs( x=expression("Stimulus Number"),
                                  y=expression('Cumulative '*Delta*'F/F'),
                                  title="",
                                  subtitle="")+
                            theme_tufte()+
                            my.theme+
                            facet_grid(eval(expr(!!ensym(secondAxis)~!!ensym(plotBy))))+
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

               			ggsave(filename='Cumulative_amplitude_vs_stimEpoch.png',plot=cumPlot, device="png",dpi=600, units="in",width=12,height=12)


						cumPlot <-ggplot(allEpochs_final, aes(x=stimEpoch,y=releaseProb)) +
														{if(all_ROIs_switch == FALSE)geom_line(aes(group=ROI_keys), alpha=0.5,colour='grey',size=1)}+
														geom_smooth(aes(group=protocol, colour=protocol, fill=protocol), formula = y~x, method="loess", se=TRUE, span=0.25)+
														#stat_summary_bin(aes(group = protocol, colour=protocol), geom="point", fun = mean, bins=numBins)+
														#stat_summary_bin(aes(group = protocol, colour=protocol), geom="line", fun = mean, bins=numBins)+
                        		#stat_summary_bin(aes(group = protocol, colour=protocol), geom="errorbar", fun = mean_se, bins=numBins)+
                        		scale_colour_brewer(palette = "Dark2")+
                            scale_fill_brewer(palette="Dark2")+
                            labs( x=expression("Stimulus Number"),
                                  y=expression('Number of Events Observed per Stimulus'),
                                  title="",
                                  subtitle="")+
                            theme_tufte()+
                            my.theme+
                            facet_grid(eval(expr(!!ensym(secondAxis)~!!ensym(plotBy))))+
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

               			ggsave(filename='ReleaseProb_vs_stimEpoch.png',plot=cumPlot, device="png",dpi=600, units="in",width=12,height=12)



} else {


	print('makeGraphs = FALSE')
}





summaryStats

}
