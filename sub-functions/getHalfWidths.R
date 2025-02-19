#getHalfWidths.R
#input fitsCleaned, output fitsCleaned with 3 columns attached.
# $rise_hw : this will be a linear interpolation between the baseline and the maximum amplitude
# $decay_hw : this will be a simple calculation of where the half-width should be based on the extracted tau 
# $half_amplitude : this is just .5 x the max amplitude. used to define display params
# within dataframe there should be a state column : $annotation. $annotation will put a different colored dot at the max and the estimated peak start index.


# in the summary dataframe, decay_hw - rise_hw will establish $halfwidth


getHalfWidths<- function(df,groupers) {


estimatedRiseTime = 0.025 #s bioRxiv says 18.9 +- 0.5 ms for GluSnFR3 peaks

interFrame = as.vector( df %>% ungroup() %>% summarise( interFrame=mean(interFrame,na.rm=TRUE) ) ) %>% 
				unlist(., use.names=FALSE)

lookBack = ceiling(estimatedRiseTime/interFrame)



	suppressMessages(getTaus<-df %>% group_by_at(groupers) %>%
					                unnest(tidied) %>%
					                spread(term,estimate) %>%
					                ungroup() %>%
					                group_by_at(groupers) %>%
					                summarise(	tau_decay = unique( tau_decay[!is.na(tau_decay)] ),																				#tau_decay
					                           A = unique( A[!is.na(A)] ) 
					                           )
					                
					    )
					

	suppressMessages(getRiseEdge<- df %>% group_by_at(groupers) %>%
					                      mutate(peakMax = case_when(dFF == max(dFF) ~ "max"),
																        peakRise = ifelse(lead(peakMax, n=lookBack) == "max", "riseEdge", NA)) %>%
					                      dplyr::filter(peakMax == "max" | peakRise == "riseEdge") %>%
					                      summarise(riseEdgeTime = absoluteTime[which(peakRise == "riseEdge")],																	#get estimated rise edge point. might be too far! 
								                     riseEdgeAmplitude = dFF[which(peakRise == "riseEdge")],
								                     maxTime = absoluteTime[which(peakMax=="max")], 																							
								                     maxAmplitude = dFF[which(peakMax=="max")])
			            )
					

	suppressMessages(halfWidths<- left_join(getRiseEdge,getTaus) %>% group_by_at(groupers) %>%
					            mutate(	     #determine rise half point
					                         rise_half_time = ( ( maxTime - riseEdgeTime )*0.5 ) + riseEdgeTime,                                                             	#linear interpolation of rise_half_time
					                         rise_half_amplitude = ( maxAmplitude )*0.5, 																						#linear interpolation of rise_half_amplitude
					                         
					                         
					                         #determine decay half point
					                         decay_half_amplitude = 0.5*(A),
					                         decay_half_time = maxTime+(-(log( decay_half_amplitude/A )*tau_decay) ),
					                         
					                         halfWidth = decay_half_time - rise_half_time,
					                         slope_of_halfpoints = decay_half_amplitude - rise_half_amplitude																#could be nice as a quality control measure? 
					                         
				                         ) 							        

					    )
					
halfWidths



}




#next steps in the pipeline
#which is better? half widths based on A or halfwidths based on true max amplitude? 
