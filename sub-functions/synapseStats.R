#def_releaseProb.R
# def_releaseProb.R is a function which will function similarly to peakStats.R
# The function will accept an input of the cleaned dataframe and output plots and statistics that summarise the activity of different ROIs with respect to their stimuli. 
# In order to properly assess the dataset, it will be necessary to accept the lookup table, stimList, as an input. 
# The desired summary statistics are as follows: 

# What was the total number of release events per ROI? (unique(peakID) ) not inclusive of "NotPeak"
# How many synchronous vs asynchronous events per ROI? compare as a % of total? Raw?
# What was the interpolated release probability per synapse? Is this worth expressing as a conditional? For example, at synapses that have a lot of async events, what was their release probability? 
# A cool graph could be dFF vs. stimulus epoch... per ROI ? 
# coefficient of variance






interSpike = function(df, stimList, groupers,groupers_byPeak, plotBy, levels, secondAxis = NULL, color_override=color_override, keys = keys, peak_keys = peak_keys,bin_override=NULL) {

								.keyvars = rlang::syms(keys)
								.peak_keyvars = rlang::syms(peak_keys)
								with_peakKeys = append(groupers_byPeak, "peakKey")
								getInterSpike<- #suppressMessages(
												df %>% mutate(stimKey = paste(!!!.keyvars, sep="-"),
																peakKey = paste(!!!.peak_keyvars, sep="-") ) %>%
														group_by_at(with_peakKeys) %>%  
													summarise(interSpike = unique(interSpike)) %>% ungroup()
												#)
								#getInterSpike<- getInterSpike[!is.na(getInterSpike$interSpike)]
								cleanInterSpike<- getInterSpike %>% ungroup() %>% group_by(timeClass) %>% mutate(interSpike_ms = interSpike*1000) %>% dplyr::filter( interSpike_ms < quantile(interSpike_ms, 0.95) )

								print("Getting dT plots")
								cleanInterSpike[[plotBy]]<- factor(cleanInterSpike[[plotBy]],
                               									levels = levels)
				


								if(is.null(secondAxis)) {
							
										print("No second axis specified.")
										var = varList
										
										source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotBoxplots.R") 
										plotBoxplots(df = cleanInterSpike, groupers = groupers ,plotBy = plotBy, var = var, color_override=color_override)

										source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotHistos.R") 
										plotHistos(df = cleanInterSpike, groupers = groupers ,plotBy = plotBy, var = var,color_override=color_override,bin_override = bin_override)

								} else {

										print("Plotting graphs by a second variable.")
										print(secondAxis)
										var = varList

										source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotBoxplots.R")
										plotBoxplots(df = cleanInterSpike, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override)

										source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/plotHistos.R") 										
										plotHistos(df = cleanInterSpike, groupers = groupers ,plotBy = plotBy, var = var, secondAxis=secondAxis, color_override=color_override,bin_override = bin_override)


									}	





							
                         cleanInterSpike
						#synapseStats

}







	#releaseProb <- suppressMessages(left_join(df,cleanInterSpike) %>% group_by_at(groupers) %>% 
								#		summarise(tau_decay_ms = unique(tau_decay)*1000,
								#					halfwidth_ms = unique(halfWidth)*1000,
								#					amplitude = max(dFF, na.rm=TRUE)) %>% na.omit() %>%   #interSpike is tricky to add because of Inf
								#			dplyr::filter(tau_decay_ms < quantile(tau_decay_ms, 0.95), 
								#							tau_decay_ms > quantile(tau_decay_ms, 0.05), 
								#							
								#							halfwidth_ms < quantile(halfwidth_ms, 0.95), 
								#							halfwidth_ms > quantile(halfwidth_ms, 0.05),  
								#							
								#							amplitude < quantile(amplitude, 0.95), 
								#							amplitude > quantile(amplitude, 0.05)#,
#
															#interSpike_s < quantile(interSpike_s, 0.95),
															#interSpike_s > quantile(interSpike_s, 0.05)
															#)
											 
										#		)


						
						#peakStats[[plotBy]]<- factor(peakStats[[plotBy]],
                         #      									levels = levels)

						







