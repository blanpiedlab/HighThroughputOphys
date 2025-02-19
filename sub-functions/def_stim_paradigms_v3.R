#def_stim_paradigms_v3.R
#First piloted on PTP-sigma dataset


#### Pseudo-code ####
# read in dataframe
# generate unique stimkey using the keyvars





def_stim_paradigms<- function(dataframe, APcount = num_APs,keys,TTL_index = TTL_start) {


					#existing protocol ids
					protocol_id1 = "BasalRP"
					protocol_id2 = "PP"
					protocol_id3 = "Hz"

					exist_protocols = c(protocol_id1,protocol_id2,protocol_id3)



					## Step 1. Generate keys for each subset of the dataframe. 
					.keyvars = rlang::syms(keys)
					stim = fitsLabeled %>% 
								mutate(stimKey = paste(!!!.keyvars, sep="-") )
					unique.stimKeys = unique(stim$stimKey)
					stimList = list() #instantiate stimList




					for (j in 1:length(unique.stimKeys)) {

								keyName = unique.stimKeys[[j]]
								j_video = stim %>% ungroup() %>% dplyr::filter(stimKey %in% keyName)
								protocol_string = unique(j_video$protocol) 									#define the protocol string to be referenced below
								
								

								#determine stimulus characteristics for each protocol

								if ( str_detect( protocol_string, protocol_id1 ) == TRUE ) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == TTL_index)])   #stimulus always start at frame 100 (so termination of 99 is a good approx) 
									
									set.sequence = first.pulse 												#stimulus is a single action potential
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a BasalRP video : video #", j, " starting at frame index = ", TTL_index))
									next

								} else if ( str_detect( protocol_string, protocol_id2 ) == TRUE ) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == TTL_index)])   
									interStimulus = as.numeric( gsub( "PP", "", as.matrix( protocol_string ) ) ) / 1000 #convert PP\\d... to \\d... and then to seconds
									second.pulse = first.pulse+interStimulus 
									
									set.sequence = c(first.pulse, second.pulse)
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a Paired Pulse video : video #", j, ", ISI was ", interStimulus, " starting at frame index = ", TTL_index))
									next

								} else if (str_detect( protocol_string, protocol_id3 ) == TRUE) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == TTL_index)])   
									interStimulus = 1 / as.numeric( gsub( "Hz", "", as.matrix( protocol_string ) ) )  #convert \\dHz... to \\d... and then take inverse to get interstimulus in seconds
									last.pulse = first.pulse+( (APcount-1) * interStimulus )
									
									set.sequence = seq(from=first.pulse, to=second.pulse, by=interStimulus)
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a stimulus train video with ", APcount, "APs : video #", j, ", train frequency was ", protocol_string, " starting at frame index = ", TTL_index))
									next

								} else {

									print( paste0( "Looks like there wasn't a protocol match in video #", j, "with keyName : ", keyName) )
									

								}

						}

					stimList

					}

