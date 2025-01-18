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

png_plotter_ROIoverlay<- function(df = df, png_dir = png_dir,prefix_str=prefix_str, suffix_str = suffix_str, dir_regex = dir_regex, marker_map = marker_map, activity_map=activity_map, save_bool=save_bool,tag_var=taggy, user_x = user_x, user_y = user_y,guide_title=guide_title,guide_bool=guide_bool,scalefactor=scalefactor){

				#fix vid_key
				fixed_vid_key = gsub("-","_",unique(df$formatted_vid_key) ) 
				#get Ca2+
				Ca_key = unique(df$Ca)
				#get protocol
				protocol_key = unique(df$protocol)

				#Derive the filename to look for. 
				
				png_dir <- png_dir
				setwd(png_dir)
				# Specify the file name you want to find
					

				# List all directories (including subdirectories) in the parent directory
				#all_dirs <- list.dirs(png_dir, full.names = TRUE)

				#Derive the .png filename to look for.
				png_filename = paste0(prefix_str,fixed_vid_key, suffix_str)  
				png_file = list.files(pattern=png_filename, full.names=FALSE, recursive=TRUE) 
				found_png = !is.null(png_file)
				print(paste0("Filename should be: ", png_filename))
				
if(found_png==TRUE){
				######read_png############ so you can change directory to read in Ch2Results.csv
								y <- readImage(png_file)
								g <- rasterGrob(y, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
								file_dir = dirname(png_file) 

								#get .csv file with ROI coordinates
								setwd(file_dir)
								csv_path = "Ch2Results.csv"				
								csv_filename<- list.files(pattern=csv_path, full.names=TRUE, recursive=TRUE) #Recurse through all subdirectories in the directory
								roi_info <- readr::read_csv(csv_filename, show_col_types = FALSE)
								colnames(roi_info)[1]<-"ROI"
								pixel_to_um = 9.9533 #pixels / 1 um
								roi_info <- roi_info %>% mutate(X_um = X*10/pixel_to_um, 
								                                Y_um = Y*10/pixel_to_um,
								                                num_padded = str_pad(ROI, 4, pad= "0"),
								                                ROINumber = paste0("ROI",num_padded))
								
								
								library(RImageJROI)
								
								ij_ROIs_marker = list.files(pattern=marker_map, full.names=TRUE, recursive=FALSE)
								ij_ROIs_activity = list.files(pattern=activity_map, full.names=TRUE, recursive=FALSE)

								ijzip_marker <- read.ijzip(ij_ROIs_marker)
								ijzip_activity <- read.ijzip(ij_ROIs_activity)
								
								

								ROIs_combined_marker = data.frame()
								ROIs_combined_activity = data.frame()
								ziprange_marker = length(ijzip_marker)
								ziprange_activity = length(ijzip_activity)

								for (i in 1:ziprange_marker) {
								
									ijROI <- ijzip_marker[[i]]
									coords_x <- rev(ijROI$coords[,1])
									coords_y <- rev(ijROI$coords[,2])
									out<- data.frame(x_pixel = coords_x, y_pixel = coords_y,ROINumber = paste0("ROI",str_pad(i, 4, pad= "0")))
									ROIs_combined_marker <- rbind(ROIs_combined_marker, out)
								}

								for (i in 1:ziprange_activity) {
								
									ijROI <- ijzip_activity[[i]]
									coords_x <- rev(ijROI$coords[,1])
									coords_y <- rev(ijROI$coords[,2])
									out<- data.frame(x_pixel = coords_x, y_pixel = coords_y,ROINumber = paste0("ROI",str_pad(i, 4, pad= "0")))
									ROIs_combined_activity <- rbind(ROIs_combined_activity, out)
								}
									
									
							  segmentation_levels = c("marker","activity")
                        
							  roi_polygons_marker<-suppressMessages(ROIs_combined_marker %>% 
							  										mutate(X_um = x_pixel/pixel_to_um, 
							                                        		Y_um = y_pixel/pixel_to_um,
							                                        		segmentation = as.character("marker")),

							  										)
							  roi_polygons_marker$segmentation <- factor(roi_polygons_marker$segmentation, levels = segmentation_levels)

						
							  roi_polygons_activity<-suppressMessages(ROIs_combined_activity %>% 
							  										mutate(X_um = x_pixel/pixel_to_um, 
							                                        		Y_um = y_pixel/pixel_to_um,
							                                        		segmentation = as.character("activity")),

							  										)
							  roi_polygons_activity$segmentation <- factor(roi_polygons_activity$segmentation, levels = segmentation_levels)


								#roi_polygons_merge = bind_rows(roi_polygons_marker,roi_polygons_activity)
						

							 scalebar = data.frame(x_start = user_x,
							 										x_end = user_x+5,
							 										y_start = user_y,
							 										y_end = user_y)

							  #axis_bool = TRUE
							  upper_range_x = 256/pixel_to_um   #ifelse(range(roi_polygons$x_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um
							  upper_range_y = 256/pixel_to_um   #ifelse(range(roi_polygons$y_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um 
				
					  	  	 source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
						




							  if(guide_bool == TRUE) {
												
							  	p1<- ggplot() +
								  annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
								  geom_polygon(data=roi_polygons_activity, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill="segmentation",colour="segmentation"),alpha=0.7,size=1)+
								  geom_polygon(data=roi_polygons_marker, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill="segmentation",colour="segmentation"),alpha=0.4,size=1)+								  
								  #geom_segment(data=scalebar, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= scalefactor*10, colour="white")+
								  
								  labs(#x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), #title = "ROI Overlay",
  											tag=tag_var)+#paste0(fixed_vid_key,": Spatial Plot with ROI overlay only")) +
								  theme_minimal()+
								  scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
								  scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+	
								  scale_colour_manual(values = c('marker' = 'blue',"activity" = "red"), guide= "none")+
								  scale_fill_manual(values = c( 'marker' = 'blue', 'activity' = "red"),labels = c('marker' = 'marker', 'activity' = 'activity'),limits = c("marker", "activity"))+
					                                            
								  
								  my.theme_spatialROI+
								   guides(fill=guide_legend(title=guide_title,override.aes = list(alpha=1, size = scalefactor*1.2)))+
								  theme(legend.justification = c(0,0), legend.position = c(1.1, 0.5),
								  		#legend.box.background = element_rect(fill = "white", color = "black"),
										axis.text = element_blank(),
										axis.title=element_blank(),
										legend.title = element_text(colour="black", size=scalefactor*28, family="sans")
										)

								  
								
								
								} else if(guide_bool == FALSE){
									p1<- ggplot() +
										annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
										geom_polygon(data=roi_polygons_activity, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill="segmentation",colour="segmentation"),alpha=0.75,size=1)+
										geom_polygon(data=roi_polygons_marker, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill="segmentation",colour="segmentation"),alpha=0.5,size=1)+
										  
										#geom_segment(data=scalebar, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 10, colour="white")+
										  
										labs(#x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), #title = "ROI Overlay",
		  										tag=tag_var)+#paste0(fixed_vid_key,": Spatial Plot with ROI overlay only")) +
										 theme_minimal()+
										 scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
										 scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+	
										 scale_colour_manual(values = c("activity" = "red", 'marker' = 'blue'), guide= "none")+
										 scale_fill_manual(values = c("activity" = "red", 'marker' = 'blue'))+
							                                           
										 
										 my.theme_spatialROI+
										 guides(colour = guide_legend(show = FALSE),
												fill = guide_legend(show = FALSE) )+
										 theme(legend.position="none",
										       axis.text = element_blank(),
												axis.title=element_blank(),
												plot.title=element_blank())
												
												}  

										  
										
										ggsave(filename=paste0(fixed_vid_key,"_",Ca_key,"Ca_spatial_plot_onlyROIs_labeled.png"),plot=p1, device="png",dpi=600, bg="white",units="in",width=16,height=16)
														# 		 	p2<- ggplot() +
												# 			  annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
												# 			  geom_polygon(data=roi_polygons, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill=var,colour=var),alpha=1,size=4)+
												# 			  scale_fill_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=guide_limits,name=guide_title, oob=scales::squish)+
												# 			  scale_colour_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=guide_limits, oob=scales::squish,guide="none")+
												# 			  geom_segment(data=ROI_ID_A, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  				# 			  geom_segment(data=ROI_ID_B, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  				# 			  geom_segment(data=scalebar, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 5, colour="white")+
								  
												# 			  #scale_color_manual(values = 'black', labels = 'Missing value') +
											  	# 			  #geom_text(data=roi_info_mutated, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",label="PPR_adj"),  colour="black", size=5, fontface='bold')+
															  
												# 			  labs(x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), title = paste0(Ca_str,", ",protocol_str),tag=tag_var  ) +  #," : ",fixed_vid_key)) +
												# 			  theme_minimal()+
												# 			  scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
												# 			  scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+
												# 			  my.theme_spatialROI+
												# 			  guides(colour = guide_legend(show = FALSE),
												# 					fill = guide_legend(show = FALSE) )+
												# 			  theme(legend.position="none",
												# 			  						axis.text = element_blank(),
												# 							  		axis.title=element_blank(),
												# 							  		plot.title=element_blank())
															  
												# 			}  
												  
												
												
												#	ggsave(filename=paste0(Ca_key,"_",protocol_key,"_",fixed_vid_key,"_",var,"_spatial_plot_.png"),plot=p2, device="png",dpi=600, bg="white",units="in",width=16,height=16)
												




	} else { print("something strange happened")}
p1
}

