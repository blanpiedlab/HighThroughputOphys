#peakStats_spont.R
#1.10.23




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

					   	onlyCleanPeaks<- suppressMessages(df %>% group_by_at(groupers) %>% dplyr::filter(peakID != "NotPeak") %>%
					   										summarise(amplitude = max(dFF))
					   									) 
					   	onlyCleanPeaks$stimEpoch[which(is.na(onlyCleanPeaks$stimEpoch))]<- 0 #coerce NAs might be necessary? 

#					   	
					   	peakStats<- suppressMessages(left_join(onlyCleanPeaks, taus) %>% group_by_at(groupers) %>% 
					   											summarise(tau_decay_ms = tau_decay*1000,
					   														amplitude = amplitude#,
					   														) %>%
					   											ungroup() %>%
					   											group_by(Ca) %>% 
					   											 dplyr::filter(	tau_decay_ms >= 10, 
					   											  				#tau_decay_ms < quantile(tau_decay_ms, 0.95, na.rm=TRUE), 
																			  	#tau_decay_ms > quantile(tau_decay_ms, 0.05, na.rm=TRUE),
																			  	tau_decay_ms <= 200,
																			  	amplitude >= 0.05,
																			  	t_half > 0
																			  	#amplitude < quantile(amplitude, 0.95, na.rm=TRUE), 
																			  	#amplitude > quantile(amplitude, 0.05, na.rm=TRUE)) 
					   									)


					   	peakStats
						
						# peakStats[[plotBy]]<- factor(peakStats[[plotBy]],
                        #        									levels = levels)

						
						 


						# if(secondAxis_switch == FALSE) {
							

						# 	print("No second axis specified.")
							

						# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotBoxplots.R")

							
						# 	for (i in 1:length(varList)) { 

						# 	var = varList[i]
						# 	plotBoxplots(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, color_override=color_override)

						# 	}	




						# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotHistos.R") 
							
							
						# 	for (j in 1:length(varList)) { 

						# 	var = varList[j]
						# 	plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var,color_override=color_override)

						# 	}	


						# } else {

						# 	print("Plotting graphs by a second variable.")
						# 	print(secondAxis)

							
						# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotBoxplots.R")

							
						# 	for (i in 1:length(varList)) { 

						# 	var = varList[i]
						# 	plotBoxplots(df = peakStats, groupers = subgroup ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

						# 	}	




						# 	source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/plotHistos.R") 
							
							
						# 	for (j in 1:length(varList)) { 

						# 	var = varList[j]
						# 	plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

						# 	}	

							

						

						# }


							
						

						
						# peakStats


}











