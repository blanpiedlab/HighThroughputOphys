## Goal: define sync vs. async 
# function structure
# accept a dataframe, grouped by groupers
# return a dataframe with sync/async definitions (timeClass) + colors (timeColor)


def_sync_async <- function(dataframe, groupers_plus,keys) {

					.keyvars = rlang::syms(keys)

					minpositive = function(x) min(x[x > 0 ])
					flexible_interSpike<- function(peakTime, stimReference) {
										stims<- stimList[[stimReference]]
										interSpike<- minpositive(peakTime - stims)
										interSpike
					}
					
					sync_async<- fitsLabeled %>% group_by_at(groupers_plus) %>%
										dplyr::filter(peakID != "NotPeak") 	%>%
										slice(which.max(dFF)) 				%>%
										mutate(timeRef = paste(!!!.keyvars, sep="-")
												)						
							
					sync_async$interSpike<- mapply(flexible_interSpike, peakTime = sync_async$absoluteTime, stimReference = sync_async$timeRef)		

					#groupers_plus_timeRef = append(groupers_plus, "timeRef")
					sync_async <- suppressMessages(sync_async %>% group_by_at(groupers_plus) %>%
									summarise(interSpike = unique(interSpike),
												timeClass = case_when(	interSpike == Inf ~ "spontaneous",
																	interSpike <= 0.030 ~ "sync",
																	interSpike > 0.030 ~ "async"
																	)
											)
                                    )

					suppressMessages(sync_async_fitsLabeled<- left_join(dataframe, sync_async)
										)
					sync_async_fitsLabeled$timeClass[is.na(sync_async_fitsLabeled$timeClass)] <- "noPeak"

					sync_async_fitsLabeled


}