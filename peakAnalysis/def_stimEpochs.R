#def_stimEpochs.R

#def_stimEpochs.R accepts a dataframe and a few string vectors containing column names and produces a new dataframe with stimulus epochs defined.

#This script should contain a few ifelse statements changing the behavior depending on which protocol dataset is being fed to the function.
#it should be noted that the column_vars "stimkey" and "ROI_keys" are being sunset in this version, as they already exist in the dataset at this point in the analysis script.


#if the unique(protocol) == "BasalRP", then we should consider there to be only one stimulus epoch.  
#if the unique(protocol) == "PP", then we should consider there to be 2 stimuli
#if 'train', then we are assuming there are 25 stimuli


def_stimEpochs = function(df = df, stimList = stimList, groupers = groupers, keysEpoch=keysEpoch) {

					.keyvars_epoch = rlang::syms(keysEpoch)
					protocol_string = unique(df$protocol)

					protocol_id_BasalRP = "singleAP"
					protocol_id_PP = "PP"
					protocol_id_TRAIN = "Hz"


					stim.df = stack(stimList)

					names(stim.df)[names(stim.df)=="ind"] <- "stimKey"
					names(stim.df)[names(stim.df)=="values"] <- "stimulusTimes"


					if (str_detect(first(protocol_string), protocol_id_PP) == TRUE 
     									& protocol_id_BasalRP %in% protocol_string) 
										{ print("There seems to be a singleAP protocol along with multiple paired pulse conditions in this dataset. Attempting a different route toward defining stimEpochs.")


										print("writing stim epoch_keys for BasalRP protocol.")
										number_of_stimuli = 1 #1 APs 


										stim.df <- suppressMessages( stim.df %>% group_by(stimKey) %>% summarise( firstStim = min(stimulusTimes),
																													lastStim = max(stimulusTimes) )
																	)

										singleAP.df <- suppressMessages( left_join(df, stim.df) %>% 
																	ungroup() %>%
																	group_by_at(groupers) %>% 
																	dplyr::filter(protocol == "singleAP") %>%
																	mutate(stimEpoch = case_when(absoluteTime >= unique(firstStim) ~ 1 ),
																			epoch_keys=paste(!!!.keyvars_epoch,sep='-'),
																			normTime = absoluteTime - firstStim )
																	)

										
										#singleAP.df

										print("writing stim epoch_keys for PairedPulse protocols.")
										#number_of_stimuli = 1 #1 APs #deprecated in paired pulse 


										
										 PP.df <- suppressMessages( left_join(df, stim.df) %>% 
																	ungroup() %>%
																	group_by_at(groupers) %>% 
																	dplyr::filter(protocol != "singleAP" ) %>%
																	mutate(stimEpoch = case_when(absoluteTime >= unique(lastStim) ~ 2,
																									absoluteTime >= unique(firstStim) ~ 1),

																			epoch_keys=paste(!!!.keyvars_epoch,sep='-'),
																			normTime = absoluteTime - firstStim )
																	)

										
										#PP.df


										all_stimKeyed.df <- rbind(PP.df,singleAP.df )
										print("The new manipulation may have worked!?")
										all_stimKeyed.df






								} else if(str_detect( first(protocol_string), protocol_id_BasalRP ) == TRUE ) {
										
										print("writing stim epoch_keys for BasalRP protocol.")
										number_of_stimuli = 1 #1 APs 


										stim.df <- suppressMessages( stim.df %>% group_by(stimKey) %>% summarise( firstStim = first(stimulusTimes) ) 
																	)

										onlyStims.df <- suppressMessages( left_join(df, stim.df) %>% 
																	ungroup() %>%
																	group_by_at(groupers) %>% 
																	dplyr::filter(absoluteTime >= unique(firstStim) ) %>%
																	mutate(stimEpoch = ntile(index, number_of_stimuli ),
																			epoch_keys=paste(!!!.keyvars_epoch,sep='-') )
																	)

										all_stimKeyed.df <- suppressMessages(left_join(df, onlyStims.df) ) 
										all_stimKeyed.df
								} else { 

									print("I guess there wasn't a suitable protocol for defining stim Epochs.")
								}
								


}
