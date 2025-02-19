#set_dataDirectory.R

checkDir<-getwd()

if (checkDir == dataDir){
  print("We are already in the desired data directory.")
  print( getwd() )
} else {
  print("Changing working directory back to the data directory.")
  setwd(file.path(dataDir) )
  print( getwd() )

}
