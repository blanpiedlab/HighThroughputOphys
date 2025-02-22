#fileRename.R
#Written by Samuel T. Barlow 

# fileRename.R carries out two functions: 
# First, it renames all .csv files in a destination as an ROI number with 4 digits. This makes file ordering during analysis more sensible. 
# Second, each .csv file is read into R as a data.frame, a few columns are appended, and the modified data.frame is then saved into a new subdirectory with the suffix "_stamped". 


rename_csv<- function(folder) {
              path=setwd(folder)
              path=setwd(folder) ##Calling this twice in a row fixes a weird bug where it lags the folder rename by 1 because it's calling the primary_folder as the initial subdirectory
              
              filenames<- list.files(pattern="ROI", full.names=TRUE,  recursive=F)
              details = file.info(filenames)
              details = details[with(details, order(as.POSIXct(mtime))),] #This line orders the files by the time they were created (ROI indexing)
              files = rownames(details)

              placeholder = basename(path)
            
              b<- sprintf("%s_ROI%04d.csv", placeholder,seq(files)) #Here, sprintf is adding padding zeros to normalize the data read
              file.rename(files,b)

  
}


read_and_write_csv<- function(file,folder, exposureTime){
                  #path=setwd(folder)
                  #path=setwd(folder)
                  #filenames<- list.files(pattern="ROI",full.names=TRUE,recursive=F)
                  folder_name_fixed<- str_remove(folder, pattern = pattern1)
                  source_info<-basename(folder)
                  source_info_fixed<-str_remove(source_info, pattern = pattern1)
                  print(paste0("This is the new foldername for saving out .csv files: ", folder_name_fixed))
          
                  newFolder <- paste0(folder_name_fixed,"_stamped/")  #create the save folder for these stamped ROIs
                  dir.create(newFolder)
                  
                  #print(paste0("This is the source info ", source_info))
                  
                  tmp_file<- read.csv(file = file,header=TRUE) 
                  colnames(tmp_file)[1]<- "Frame"
                  colnames(tmp_file)[2]<- "Intensity"
                  tmp_file$fileID<- source_info_fixed
                  tmp_file$index <-seq.int(nrow(tmp_file))
                  tmp_file$timeStamps<- "placeholder"
                  tmp_file$dt <- tmp_file$index*exposureTime - exposureTime
                  #print(head(tmp_file))
                  
                  newfilename<- paste0(newFolder,str_extract(file, pattern=newFilePat),"stamped.csv")
                  
                  write.csv(x = tmp_file, file = newfilename,row.names=FALSE)
                                           #read in the stampfile that we've kept
                  
}





primary_folder = dir
subfolders = grep(pattern="_Results",list.dirs(primary_folder,  full.names = T, recursive = T), value = T)


lapply(subfolders, rename_csv)



for (i in subfolders){
      tmp_subfolder<- i
      setwd(tmp_subfolder)
      setwd(tmp_subfolder)
      fnames <- list.files(pattern="ROI", full.names=TRUE,recursive=FALSE)      
      
      lapply(fnames, read_and_write_csv,folder = tmp_subfolder, exposureTime = frameTime)
      
}

    
