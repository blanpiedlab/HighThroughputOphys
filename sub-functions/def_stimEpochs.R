#def_releaseProb.R
# def_releaseProb.R is a function which will function similarly to peakStats.R
# The function will accept an input of the cleaned dataframe and output plots and statistics that summarise the activity of different ROIs with respect to their stimuli. 
# In order to properly assess the dataset, it will be necessary to accept the lookup table, stimList, as an input. 
# The desired summary statistics are as follows: 



# What was the total number of release events per ROI? (unique(peakID) ) not inclusive of "NotPeak"
# How many synchronous vs asynchronous events per ROI? compare as a % of total? Raw?
# To express these observations (n_tot, n_sync, n_async) as a release probability, it will be necessary to define stimulus epochs 
# This should be possible using a similar approach to peakFinder_Lite.R:
# 


# What was the interpolated release probability per synapse? Is this worth expressing as a conditional? For example, at synapses that have a lot of async events, what was their release probability? 
# A cool graph could be dFF vs. stimulus epoch... per ROI ? 
# coefficient of variance



#Pseudo-code for defining stimulus epochs
# code should be able to check: 1. how many stimuli there are in an arbitrary stimList[[stimKey]]. 
# 								2.  





def_stimEpochs = function(df, stimList,groupers, 
									keys=keys, 
									keysROI=keysROI,
									keysEpoch=keysEpoch
									) {

								.keyvars = rlang::syms(keys)
								.keyvars_ROI = rlang::syms(keysROI)
								.keyvars_epoch = rlang::syms(keysEpoch)
					

								number_of_stimuli = 25 #25 APs 


								stim.df = stack(stimList)
								names(stim.df)[names(stim.df)=="ind"] <- "stimKey"
								names(stim.df)[names(stim.df)=="values"] <- "stimulusTimes"
								stim.df <- suppressMessages( stim.df %>% group_by(stimKey) %>% 
																			summarise( firstStim = first(stimulusTimes),
																						lastStim = last(stimulusTimes) 
																						) 
																		)


								stimKeyed.df <- df %>% ungroup() %>% mutate(stimKey = paste(!!!.keyvars, sep="-"))
								onlyStims.df <- suppressMessages( left_join(stimKeyed.df, stim.df) %>% 
															ungroup() %>%
															group_by_at(groupers) %>% 
															dplyr::filter(absoluteTime >= unique(firstStim),
																			absoluteTime <= unique(lastStim)
																							) %>%
															mutate(stimEpoch = ntile(index, 
																					number_of_stimuli ),
																	ROI_keys = paste(!!!.keyvars_ROI, sep="-"),
																	epoch_keys=paste(!!!.keyvars_epoch,sep='-')
																	)   																					# a bit of a tricky expression. Here we are finding the number of stimulusEpochs within each group. 
																																																		  # First, we extract the start time and define it per group. Then we extract only the section of the trace where stimuli are occurring:
																																																		  # i.e. from the start time (firstStim) to the final time (firstStim + 25 seconds)
																																																		  # Once we have only the section of each trace which corresponds to stimuli, we can n_rank (ntile()) our absoluteTime.
																																																		  # To determine the proper number of n_ranks per group (the number of stimulus epochs), we take the difference between our first time and our last time and round up - this may produce unexpected behavior if the stimulus epoch failed to be exactly 25 seconds long
																																																		  # Then we can determine the correct number of stimuli by dividing by the inverse of the stim frequency - e.g. 2 Hz for 25 seconds is 50 stimuli .... 1/2 stimuli / second = 0.5 seconds / stimulus, 25 seconds/ 0.5 seconds/ stimulus = 50 stimuli 
																					
																	)

															

								
								all_stimKeyed.df <- suppressMessages(left_join(df, onlyStims.df) ) 
								all_stimKeyed.df
								

}
