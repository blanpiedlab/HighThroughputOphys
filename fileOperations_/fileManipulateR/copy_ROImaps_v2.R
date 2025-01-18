
#copy_ROImaps_v2.R
#Written by Samuel T. Barlow
# 3/27/23

#copy_ROImaps.R is a function that will accept 4 inputs: 
# zstack_dir: a directory in which to find ROImap.zip for copying. zstack_dir is necessary to find the files to be copied, and determine which folders it can be copied to. 
# zstack_pat: a string matching pattern by which to find the ROImap.zip
# vid_dir: a directory in which to find all the .tif files that will serve as a reference for where to copy ROImap.zip
# vid_pat: a string matching pattern by which to identify the correct files. 

#copy_ROImaps.R will take these 4 inputs and produce one output: 
# The successful distribution of the correct ROImap.zip to the correct video files' folders, such that the ~20 videos or so that will have the same ROImap receive the correct file. 
# we will do this in the following order: 
### read in the list of possible files with their full.names.
### in the case of ROImaps, extract a file pattern to match the beginning of each video name
### in the case of vid_files, extract the video files and figure out their directories

### intersect the lists based on the file pattern to match the beginning of each video name
### for each ROImap.zip, find vidfiles that return true on str_detect, determine their sub-directories (subDirs)
### file.copy(pathFrom=ROImap.zip[i], pathTo=subDirs )

allpatterns <- function(fnames, patterns) {
  		i <- sapply(fnames, function(fn) all(sapply(patterns, grepl, fn)) )
  		fnames[i]
  		
}


copy_ROImaps<- function(zstack_dir = stack_dir, zstack_pat = stack_pat, video_dir = vid_dir, video_pat=vid_pats, matches = match_pats,dirname_sort = TRUE) {



					ROImaps = list.files(path=zstack_dir, pattern=zstack_pat, full.names=TRUE, recursive=TRUE)
					
					
					vid_files = list.files(path=video_dir, pattern=".tif$", full.names=TRUE, recursive=TRUE)
					filenames_out = allpatterns(vid_files, video_pat)

					for (ROImap in ROImaps) {
							if(dirname_sort == FALSE) {#if dirname_sort is false, we are not able to use the folder name to do ROImap matching. Instead, find the names of sister files in the ROImap.
								current_dir = dirname(ROImap)
								sister_vid = list.files(path=current_dir, pattern="GluSnFR3_(.*?).tif$", full.names=TRUE, recursive=TRUE)
								match_pat = str_extract(sister_vid, matches) 
							} else {
								current_dir = dirname(ROImap)
								match_pat = str_extract(ROImap, matches)
							}
							
							

							keeps = filenames_out[ which( str_detect(filenames_out, match_pat) == TRUE )  ] 
							to_subDirs = dirname(keeps)

							path_in = paste0(current_dir, '/', basename(ROImap) ) ### current method, unclear if this will work for all cases but seems promising
							#path_in = ROImaps[which(str_detect(ROImaps, match_pat) == TRUE)]
							path_out = paste0(to_subDirs, '/', basename(path_in))
							file.copy(from=path_in, to=path_out, copy.mode = TRUE, copy.date = FALSE)   
							  
							
					}


}
