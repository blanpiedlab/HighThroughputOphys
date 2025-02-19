#write_output_csv.R

#Written by Samuel T Barlow
#12.8.23


#write_output_csv.R accepts a list of dataframes. FOr each element in the list, it writes a .csv file with the list item's name. 


subDir <- "csv_output"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 




write_output_csv<- function(list_df = output_data, prefix = NULL){

					if(!is.null(prefix)){
						title_switch = TRUE
					} 


					for (i in 1:length(list_df)){
						    tmp_df <- list_df[i]
						    if(title_switch){
						    	csv_title = paste0(prefix,names(tmp_df),".csv")
						    }else{csv_title = paste0("final_output_",names(tmp_df),".csv")}


						    write.csv(tmp_df, file = csv_title)
						}
					}