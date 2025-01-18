#def_stim_paradigms_v4.R
#Written by Samuel T. Barlow 10.25.22

# Changed the stimList calculation to account for modified protocol
# all TTLs are delayed by 50 ms upon arriving at the Grass S88X
# Also I think that TTL_index - interframe - 10 ms will better approximate the first TTL


#### Pseudo-code ####
# read in dataframe
# generate unique stimkey using the keyvars





def_stim_paradigms<- function(dataframe, APcount = num_APs,keys,offset=NULL) {

					dur_TTL = 0.01

					if(!is.null(offset)){
						time_shift = offset

					} else {
						time_shift = 0
					}
					#existing protocol ids
					protocol_id_BasalRP = "singleAP"
					protocol_id_PP = "PP"
					protocol_id_TRAIN = "Hz"

					exist_protocols = c(protocol_id_BasalRP,protocol_id_PP,protocol_id_TRAIN)



					## Step 1. Generate keys for each subset of the dataframe. 
					.keyvars = rlang::syms(keys)
					stim = dataframe %>% 
								mutate(stimKey = paste(!!!.keyvars, sep="-") )
					unique.stimKeys = unique(stim$stimKey)
					stimList = list() #instantiate stimList




					for (j in 1:length(unique.stimKeys)) {

								keyName = unique.stimKeys[[j]]
								j_video = stim %>% ungroup() %>% dplyr::filter(stimKey %in% keyName)
								protocol_string = unique(j_video$protocol) 									#define the protocol string to be referenced below

								usual_frame_time<-median(j_video$interFrame,na.rm=TRUE)
								frame_jitter = dur_TTL+usual_frame_time
								

								#determine stimulus characteristics for each protocol

								if ( str_detect( protocol_string, protocol_id_BasalRP ) == TRUE ) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == (unique(j_video$TTL_start)-1) )]) + time_shift #- frame_jitter # 10 ms   
									
									set.sequence = first.pulse 												#stimulus is a single action potential
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a BasalRP video : video #", j, " starting at frame index = ", unique(j_video$TTL_start)))
									next

								} else if ( str_detect( protocol_string, protocol_id_PP ) == TRUE ) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == (unique(j_video$TTL_start)-1) )]) + time_shift #- frame_jitter
									interStimulus = as.numeric( gsub( "PP", "", as.matrix( protocol_string ) ) ) / 1000 #convert PP\\d... to \\d... and then to seconds
									second.pulse = first.pulse+interStimulus 
									
									set.sequence = c(first.pulse, second.pulse)
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a Paired Pulse video : video #", j, ", ISI was ", interStimulus, " starting at frame index = ", unique(j_video$TTL_start)))
									next

								} else if (str_detect( protocol_string, protocol_id_TRAIN ) == TRUE) {

									first.pulse = unique(j_video$absoluteTime[which(j_video$index == (unique(j_video$TTL_start)-1) )]) + time_shift #- frame_jitter
									interStimulus = 1 / as.numeric( gsub( "Hz", "", as.matrix( protocol_string ) ) )  #convert \\dHz... to \\d... and then take inverse to get interstimulus in seconds
									last.pulse = first.pulse+( (APcount-1) * interStimulus )
									
									set.sequence = seq(from=first.pulse, to=last.pulse, by=interStimulus)
									stimList[[keyName]] <- set.sequence

									print(paste0("Entered a stimList entry for a stimulus train video with ", APcount, "APs : video #", j, ", train frequency was ", protocol_string, " starting at frame index = ", unique(j_video$TTL_start)))
									next

								} else {

									print( paste0( "Looks like there wasn't a protocol match in video #", j, "with keyName : ", keyName) )
									

								}

						}

					stimList

					}

