#subtract_by_index.R

#written 10/9/23 STB

#The goal of this script is to average traces together. It accepts a dataframe containing all data as an input, along with a vector containing the grouping variables to be used. 

#Then, it re-aligns all the data via a new index (defined by which.min() to find the index closest to normTime = 0).

#finally, it recalculates a PPR_data that is output as a dataframe to be fed into a plotting function in PPR_plotter_v1.R 




subtract_by_index<- function(df=df,groupers=groupers,groupers_avg=groupers_avg, ntile_bins = bins_for_ntile){


				#reindex the dataset by the index closest to time 0
 				#df_reindex <- suppressMessages(left_join(df, df_timezero_idx)  %>% group_by_at(groupers) %>% mutate(reindex = index - timezero_idx,
 				#																									new_normTime = absoluteTime - absoluteTime[which(reindex == 0)],
 				#																									ntile_time = ntile(new_normTime, n=ntile_bins) ) ) 

				#define which index is closest to normalized time 0
				df_timezero_idx<- suppressMessages(df %>% group_by_at(groupers) %>% summarise(timezero_idx = index[which.min(abs(0-normTime))]) )

				#reindex the dataset by the index closest to time 0
 				df_reindex <- suppressMessages(left_join(df, df_timezero_idx)  %>% group_by_at(groupers) %>% mutate(reindex = index - timezero_idx) %>% arrange(reindex) )

 				#group by "reindex" so we can calculate an average of all the replicates
 				groupers_avg_revised = c(groupers_avg,"reindex")

 				#calculate the average trace
 				df_avg_by_repl <- suppressMessages(df_reindex %>% group_by_at(groupers_avg_revised) %>% summarise(avg_normTime = mean(normTime),avg_dFF = mean(dFF)) %>% arrange(reindex)) 
 				

 				#get just the singleAP trials
 				singleAP_avg <- df_avg_by_repl %>% dplyr::filter(protocol == "singleAP",reindex >= 0) 

 				#get just the columns we need to merge with PP dataset
 				singleAP_avg_kept = subset(singleAP_avg, select = c(Ca,trackROI_key,reindex,avg_dFF)) 

 				#rename the avg_dFF column to singleAP_dFF
 				singleAP_avg_v0 = singleAP_avg_kept %>% rename("singleAP_dFF" = "avg_dFF" ) %>% group_by(Ca,trackROI_key) %>% arrange(reindex)

 				#df_subtracted = suppressMessages(left_join(df_avg_by_repl, singleAP_avg_v0) %>% group_by_at(groupers_avg) %>% mutate(dFF_subtracted = avg_dFF - singleAP_dFF) %>% arrange(reindex))
 				df_subtracted = suppressMessages(left_join(df_reindex, singleAP_avg_v0) %>% group_by_at(groupers) %>% mutate(dFF_subtracted = dFF - singleAP_dFF) %>% arrange(reindex))
 				
 				#df_subtracted = df_subtracted %>% select(Ca,protocol,trackROI_key,reindex,avg_normTime,avg_dFF,dFF_subtracted)
 				df_subtracted = unique(df_subtracted)
 				#df_subtracted = suppressMessages(left_join(df_reindex,df_subtracted) %>% group_by(Ca,protocol,vid_key,trackROI_key,ROI_key,exposeNum) %>% arrange(reindex))

#df_timezero_idx
#df_reindex
 #df_avg_by_repl
 #singleAP_avg_v0
df_subtracted
}
                                                                        