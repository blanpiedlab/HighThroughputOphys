#peakStats.R




peakStats<- function(df,dfHalfWidths,groupers,plotBy, varList, levels, secondAxis = NULL, color_override=color_override){
						
						secondAxis_switch = !is.null(secondAxis)
						.group_vars = ifelse(secondAxis_switch == TRUE, c(plotBy, secondAxis), plotBy)

						peakStats <- suppressMessages(left_join(df,dfHalfWidths) %>% group_by_at(groupers) %>% 
									summarise(tau_decay_ms = unique(tau_decay)*1000,
												halfwidth_ms = unique(halfWidth)*1000,
												amplitude = max(dFF, na.rm=TRUE)) %>% na.omit() %>%   
												ungroup() %>% 
												group_by_at(.group_vars) %>% 
											dplyr::filter(tau_decay_ms < quantile(tau_decay_ms, 0.95), 
															tau_decay_ms > quantile(tau_decay_ms, 0.05),
															tau_decay_ms > 15, 
															
															halfwidth_ms < quantile(halfwidth_ms, 0.95), 
															halfwidth_ms > quantile(halfwidth_ms, 0.05),  
															halfwidth_ms > 15,

															amplitude < quantile(amplitude, 0.95), 
															amplitude > quantile(amplitude, 0.05)#,

															#interSpike_s < quantile(interSpike_s, 0.95),
															#interSpike_s > quantile(interSpike_s, 0.05)
															)
											 
												)


						
						peakStats[[plotBy]]<- factor(peakStats[[plotBy]],
                               									levels = levels)

						


						

						if(missing(secondAxis)) {
							

							print("No second axis specified.")
							

							source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotBoxplots.R")

							
							for (i in 1:length(varList)) { 

							var = varList[i]
							plotBoxplots(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, color_override=color_override)

							}	




							source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotHistos.R") 
							
							
							for (j in 1:length(varList)) { 

							var = varList[j]
							plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var,color_override=color_override)

							}	


						} else {

							print("Plotting graphs by a second variable.")
							print(secondAxis)

							
							source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotBoxplots.R")

							
							for (i in 1:length(varList)) { 

							var = varList[i]
							plotBoxplots(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

							}	




							source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotHistos.R") 
							
							
							for (j in 1:length(varList)) { 

							var = varList[j]
							plotHistos(df = peakStats, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

							}	

							

						

						}


							
						

						
						peakStats


}











