	#png_plotter.R

#written 10/10/23 STB

#The goal of this script is to generate raster plots of neuron images from .png files and overlay statistics on top of them. 
#to do this, the script is broken down into three parts.

#in part 1, a dataframe that has already been parsed to a specific video will be input into the function. We will tell the function for which data it should be looking (e.g. vid_key) and in which directory it should be looking. 
#Through this sequence, the function will idenpngy 1) a .png file which is the max-intensity projection of the z-stack; and 2) a .csv file which is the coordinate file. 

#in part 2, the data has been idenpngied and is read into the R environment. We will generate a raster plot of just the .png and mark the ROIs with red x's and ROI labels as a data control.

#in part 3, one or more arguments connoting the data to be plotted will be input as variables into the function. This should read as an if else statement:
# if (length(argin) == 1 ) { plot just the one plot} else if (length(argin) > 1 ) { for argument in argin ( plot a plot, then move onto the next in the list)}

#With this construction, the hope is that the function can be called within a larger set of plotting/analysis functions. If it sufficiently generalized, it can be deployed across a variety of datasets. 
# in principle, all it needs is a data_directory in which to find the .png and .csv of coordinates, a list of variables, and the dataset to be measured and expressed on the raster.



#### just for plotting ROI_overlays to reduce unnecessary looping

library(ggplot2)
library(grid)
library(dplyr)
library(stringr)
library(imager)
library(readr)
library(bioimagetools)
library(EBImage)
library(ggrepel)

source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme_spatialROI.R") 

png_plotter_ROIoverlay<- function(df, png_dir,prefix_str, suffix_str, dir_regex,save_bool,tag_var, ROI_IDs){

				#fix vid_key
				fixed_vid_key = gsub("-","_",unique(df$vid_key) ) 
				getROIs = df %>% dplyr::filter(ROINumber %in% ROI_IDs)
				ROIs_to_display = unique(getROIs$ROINumber)
				title_info = paste(sort(unique(getROIs$chemical_condition)), collapse="&")
				
				print(fixed_vid_key)
				#get Ca2+
				#Ca_key = unique(df$Ca)
				#get protocol
				#protocol_key = unique(df$protocol)

				#Derive the filename to look for. 
				
				png_dir <- png_dir
				setwd(png_dir)
				print(getwd() )
				# Specify the file name you want to find
					

				# List all directories (including subdirectories) in the parent directory
				#all_dirs <- list.dirs(png_dir, full.names = TRUE)

				#Derive the .png filename to look for.
				png_filename = paste0(prefix_str,fixed_vid_key, suffix_str)  
				print(png_filename)
				png_file = list.files(pattern=png_filename, full.names=FALSE, recursive=TRUE) 
				print(png_file)
				found_png = !is.null(png_file)
				
if(found_png==TRUE){
				######read_png############ so you can change directory to read in Ch2Results.csv
								y <- readImage(png_file)
								g <- rasterGrob(y, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
								file_dir = dirname(png_file) 

								#get .csv file with ROI coordinates
								#setwd(file_dir)
								#csv_path = "Ch2Results.csv"				
								#csv_filename<- list.files(pattern=csv_path, full.names=TRUE, recursive=TRUE) #Recurse through all subdirectories in the directory
								#roi_info <- readr::read_csv(csv_filename, show_col_types = FALSE)
								#colnames(roi_info)[1]<-"ROI"
								pixel_to_um = 9.9533/2 #pixels / 1 um
								#roi_info <- roi_info %>% mutate(X_um = X*10/pixel_to_um, 
								#                                Y_um = Y*10/pixel_to_um,
								#                                num_padded = str_pad(ROI, 4, pad= "0"),
								#                                ROINumber = paste0("ROI",num_padded))
								
								
								library(RImageJROI)
								#library(spatstat)
								#library(maptools)
								#library(sf)
								
								ij_ROIs = list.files(pattern="ROImap.zip", full.names=TRUE, recursive=FALSE)

								#file <- file.path(system.file(package = "RImageJROI"), "extdata", "ijroi", "polygon.roi")
								ijzip <- read.ijzip(ij_ROIs)
								#plot(x, col = "red")

								

								ROIs_combined = data.frame()
								# #owinList <- SpatialPolygons(list())
								ziprange = length(ijzip)
								for (i in 1:ziprange) {
								# 
									ijROI <- ijzip[[i]]
									coords_x <- rev(ijROI$coords[,1])
									coords_y <- rev(ijROI$coords[,2])
									out<- data.frame(x_pixel = coords_x, y_pixel = coords_y,ROINumber = paste0("ROI",str_pad(i, 4, pad= "0")))
									ROIs_combined <- rbind(ROIs_combined, out)
								}
									
							  
							  roi_polygons<-suppressMessages(left_join(ROIs_combined,df) %>% 
							  									mutate(X_um = x_pixel/pixel_to_um, 
							                                    		Y_um = y_pixel/pixel_to_um) %>%
							  									dplyr::filter(ROINumber %in% ROIs_to_display)
							  										)

							  


							  #roi_info<- suppressMessages(left_join(roi_info,df)  %>% mutate(ROI_integer = as.numeric( gsub("ROI","",ROINumber ) ) ) 
							  									#) 																	#%>% mutate(across(all_of(vars_to_plot), signif(.x,digits=3) ) )

							  ###### arrow math stuff ####
							 				# ROI_ID_A<- roi_info %>% dplyr::filter(ROINumber == ROI_IDs[1]) %>% mutate(head_end_x = X_um,
							 				# 																			head_end_y = Y_um - 0.2,
							 				# 																			tail_start_x = X_um,
							 				# 																			tail_start_y = Y_um -2)
							 				#ROI_ID_B<- roi_info %>% dplyr::filter(ROINumber == ROI_IDs[2]) %>% mutate(head_end_x = X_um,
							 				#																			head_end_y = Y_um - 0.2,
							 				#																			tail_start_x = X_um,
							 				#																			tail_start_y = Y_um -2)

							 scalebar = data.frame(x_start =  2,
							 										x_end = 22,
							 										y_start = 49,
							 										y_end = 49)


							  #axis_bool = TRUE
							  upper_range_x = 256/pixel_to_um#ifelse(range(roi_polygons$x_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um
							  upper_range_y = 256/pixel_to_um#ifelse(range(roi_polygons$y_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um 
							  	 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
						
							  	p1<- ggplot() +
								  annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
								  geom_polygon(data=roi_polygons, aes(x=X_um, y=Y_um,group=ROINumber),colour="white",fill="white",alpha=0.7,size=2)+
								  #geom_segment(data=ROI_ID_A, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  #geom_segment(data=ROI_ID_B, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  #geom_label(data=ROI_ID_A, aes(x=tail_start_x,y=tail_start_y,label=ROINumber), colour="black",size=3,fontface="bold")+
								  #geom_label(data=ROI_ID_B, aes(x=tail_start_x,y=tail_start_y,label=ROINumber), colour="black",size=3,fontface="bold")+
								  #scale_fill_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=c(0,3),alpha=1)+
								  #scale_color_manual(values = 'black', labels = 'Missing value') +
				  				  #geom_text(data=roi_info_mutated, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",label="PPR_adj"),  colour="black", size=4, fontface='bold')+
								  geom_segment(data=scalebar, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 8, colour="white")+
								  
								  labs(#x = expression("X ("*mu*"m)"), 
								  		#	y = expression("Y ("*mu*"m)"), 
								  			#title = paste0(fixed_vid_key,title_info),
								  			tag=tag_var)+#paste0(fixed_vid_key,": Spatial Plot with ROI overlay only")) +
								  theme_minimal()+
								  scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
								  scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+	
								  #if(axis_bool){easy_remove_x_axis(what=c("title","text") )+
								  #easy_remove_y_axis(what=c("title","text") )}+
								  
								  my.theme_spatialROI


								  if(save_bool == TRUE) {
									print("this imaging region should be saved")
									

									ggsave(filename=paste0(fixed_vid_key,"_spatial_plot_onlyROIs_COLOR.png"),plot=p1, device="png",dpi=600, bg="white",units="in",width=16,height=16)
									p1
								}
				
								
								
									#ggsave(filename=paste0(fixed_vid_key,"_spatial_plot_onlyROIs.png"),plot=p1, device="png",dpi=600, bg="white",units="in",width=16,height=16)
								
								# p2<- ggplot() +
								#   annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
								#   geom_polygon(data=roi_polygons, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber"),colour="black",fill="red",alpha=0.5,size=0.5,name="ROI")+
								#   #scale_fill_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=c(0,3),alpha=1)+
								#   #scale_color_manual(values = 'black', labels = 'Missing value') +
				  				#   #geom_label(data=roi_info, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",label= "ROI_integer") , fill="white", colour="black", size=6, fontface='bold')+
								  
								#   labs(x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), title = "")+#paste0(fixed_vid_key,": Spatial Plot with ROI overlay only")) +
								#   theme_minimal()+
								#   scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
								#   scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+
								  
								#   my.theme_spatialROI
								
								
									#ggsave(filename=paste0(fixed_vid_key,"_spatial_plot_onlyROIs_labeled.png"),plot=p2, device="png",dpi=600, bg="white",units="in",width=16,height=16)
						





	} else { print("something strange happened")}

}
