#Use this script to merge all ROI files in a directory with their corresponding timestamps extracted from .ims files using recursive_ims_timeStamps.m in MatLab
#Approach uses stringr (str_detect()) to extract only the timestamp file which corresponds to the matching string [i] in toMatch.
#Currently, this script requires the user to define the strings that need to be matched.
#Future iterations of this code will retain information about the source file in the ROI file title?

  
tic()
merger<- function(file,timeStamps) 
{
  tempfile<- read.csv(file = file, header=FALSE)#,fileEncoding = "UTF-8")
  stamped<- cbind(tempfile,timeStamps)
  colnames(stamped)[1]<-"Frame"
  colnames(stamped)[2]<-"Intensity"
  if(multichannel == TRUE) {
    stamped$sensor <- sensor
  }
  
  newfilename<- paste0(newFolder,str_extract(file, pattern=newFilePat),"stamped.csv")
  write.csv(x = stamped, file = newfilename,row.names=FALSE)
  
}

primary_folder = dir
setwd(primary_folder)
toMatch <- match
for (i in toMatch) 
{                 
  
  stampList <- list.files(pattern = "_timestamps", full.names=TRUE, recursive=TRUE) #get all the timestamp files
  stampLogical<-str_detect(stampList, fixed(i))                                     #find timestamp files that have a match with this loop iteration's string in "toMatch"
  if(length(stampLogical[stampLogical == TRUE]) > 1){
    #print("EXTRA TRUES found, need to find exact match")
    stamp <- stampList[stampLogical] 
    stamp<- stamp[which(!str_detect(stamp,'Glu_\\d_'))]
  } else {
    #print("Stamp should work?")
    stamp<- stampList[stampLogical]                                                   #only keep the timestamp that gives TRUE
  }
  newFolder <- paste0(file_path_sans_ext(primary_folder),"\\",i,"_stamped\\")  #create the save folder for these stamped ROIs
  dir.create(newFolder)                                                                   

  stampfile<- read.csv(file = stamp,header=TRUE)                                          #read in the stampfile that we've kept
  stampfile<- stampfile[1:(nrow(stampfile)),]                                           #select only the rows which correspond to our sliced videos (nrow-1 for videos with coord_flip)
  allROIs<- list.files(pattern="_Results_ROI\\d\\d\\d\\d",full.names=TRUE,recursive=TRUE) #read in every ROI file (had to do this because the file patterns have regex symbols)
  ROILogical<-str_detect(allROIs, fixed(i))                                               #apply the same strategy as above to only keep ROI files that match our input, i
  ROIList<- allROIs[ROILogical]                                                           #only keep files where match to i = TRUE
  
  lapply(ROIList, merger,timeStamps = stampfile)                                          #this should create all timestamps
  

}
toc()
print("Timestamp files have now been merged with their respective intensity-time traces. Woohoo!")