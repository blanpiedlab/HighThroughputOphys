## Read files into R 



#Define read_csv_filename, which names existing columns and adds new columns according to extracted strings
csv_reader <- function(fname, prefix){
		  
		  ret <- read.csv(fname)
		  colnames(ret)[1]<- "rowNumber"
		  #output_csv_list results in paste(df_name, colnames(df)). this gsub call reverses that with a user-defined prefix. 
		  colnames(ret) <- gsub(prefix,"",colnames(ret))
		  

		  ret 		  
		}



