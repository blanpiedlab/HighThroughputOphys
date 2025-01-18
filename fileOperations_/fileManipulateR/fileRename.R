#Use this script to rename files
#Code seems to need to run twice in order to properly work.  First run gives an empty list in subfolders. Unclear why.

tic()
primary_folder = dir
subfolders = grep(pattern="_Results",list.dirs(primary_folder,  full.names = T, recursive = T), value = T)

      renameFunc<- function(folder) {
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
      
lapply(subfolders,renameFunc)
toc()

print("ROI files from ImageJ have been renamed to include padding zeros in the filename.  Now we can properly sort them! Woohoo!")