#adapted from def_sync_async.R
#Written Samuel T Barlow 12.6.23

#def_interSpike accepts a dataframe and stimList and returns a dataframe with interspike features for individual peakIDs identified. It requires a suitable grouping variable. 

def_interSpike <- function(df, groupers,keys, list_stims = list_stimuli) {

					.keyvars = rlang::syms(keys)
					
					minpositive = function(x) min(x[x > 0 ])
					flexible_interSpike<- function(peakTime, stimReference) {
										stims<- stimList[[stimReference]]
										interSpike<- minpositive(peakTime - stims)
										interSpike
					}
					
					interSpike_df<- df %>% group_by_at(groupers) %>%
										dplyr::filter(peakID != "NotPeak") 	%>%
										slice(which.max(dFF)) 				%>%
										mutate(stimKey = paste(!!!.keyvars, sep="-"),
												)						
							
					interSpike_df$interSpike<- mapply(flexible_interSpike, peakTime = interSpike_df$absoluteTime, stimReference = interSpike_df$stimKey)		

					interSpike_df <- suppressMessages(interSpike_df %>% group_by_at(groupers) %>%
									summarise(interSpike = unique(interSpike))
									)

					labeled_df<- suppressMessages(left_join(df, interSpike_df) %>% mutate(stimKey = paste(!!!.keyvars, sep="-"))
										)

					
					labeled_df


}