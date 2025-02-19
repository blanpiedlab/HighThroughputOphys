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





# # Load required packages
# library(imager)
# library(readr)

# # Define png_reader function
# png_reader <- function(directory_path, trackROI_keys) {
#   # Initialize empty lists to store data frames
#   merged_data_list <- list()
#   png_values_list <- list()
  
#   # Iterate through the keys in the trackROI_keys dictionary
#   for (key in names(trackROI_keys)) {
#     # Get file paths for the current .png and .csv files
#     png_file_path <- file.path(directory_path, paste0(trackROI_keys[[key]], ".png"))
#     csv_file_path <- file.path(directory_path, paste0(trackROI_keys[[key]], ".csv"))
    
#     # Read the .png file
#     image_data <- imager::load.image(png_file_path)
    
#     # Read ROI coordinates from the .csv file
#     roi_coordinates <- readr::read_csv(csv_file_path)
    
#     # Merge ROI coordinates with image data
#     merged_data <- merge(as.data.frame(image_data), roi_coordinates, by.x = "common_column_in_image_data", by.y = "common_column_in_roi_coordinates")
    
#     # Add merged data frame to the list
#     merged_data_list[[key]] <- merged_data
    
#     # Add .png values data frame to the list
#     png_values <- as.data.frame(image_data)
#     png_values_list[[key]] <- png_values
#   }
  
#   # Return a list containing merged data frames and .png values data frames
#   result_list <- list(merged_data = merged_data_list, png_values = png_values_list)
#   return(result_list)
# }

# # Example usage of png_reader function
# # Specify the directory path and trackROI_keys dictionary
# directory_path <- "path/to/your/directory"
# trackROI_keys <- list("File1" = "file1", "File2" = "file2")

# # Call the png_reader function
# result <- png_reader(directory_path, trackROI_keys)

# # Now, result contains a list with two elements:
# # - result$merged_data: A list of merged data frames with xy coordinates for each .png file
# # - result$png_values: A list of data frames containing .


library(ggplot2)
library(grid)
library(dplyr)
library(stringr)
library(imager)
library(readr)
library(bioimagetools)
library(EBImage)

source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme_spatialROI.R") 



# check<- c("dish05-plate10-region01",
# 			"dish05-plate10-region02",
# 			"dish05-plate11-region01",
# 			"dish05-plate11-region02",
# 			"dish06-plate01-region01",
# 			"dish06-plate01-region02",
# 			"dish06-plate02-region01",
# 			"dish06-plate02-region02") 
# png_dir = "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/imageFiles/zstack/lowpwr_zstack/zstackOutputs_SQ"
# prefix_str = "^MAX_GluSnFR3_SyPhy_"
# suffix_str = "_0pt5Ca_zstack_lowpwr__FLATTENED.png"
# dir_regex = "GluSnFR3(.*?)region\\d\\d_"
png_plotter<- function(df = df, png_dir = png_dir,prefix_str=prefix_str, suffix_str = suffix_str, dir_regex = dir_regex, 
						vars_to_plot=vars_to_plot,save_bool = save_bool,guide_bool=guide_bool,
						tag_var=tmp_tag_var, ROI_IDs = ROInames,guide_title=guide_title,guide_limits=guide_limits){



				vector.is.empty <- function(x) return(length(x) ==0)  



				#fix vid_key
				fixed_vid_key = gsub("-","_",unique(df$vid_key) ) 
				print(paste0("vid_key is : ", fixed_vid_key) )
				#get Ca2+
				Ca_key = unique(df$fix_Ca)
				print(paste0("Ca_key is : ", Ca_key) )
				#get protocol
				protocol_key = "singleAP" #unique(df$protocol)
				print(paste0("protocol_key is : ", protocol_key) )

				#Derive the filename to look for. 
				
				png_dir <- png_dir
				setwd(png_dir)
				setwd(png_dir)
				
				#Derive the .png filename to look for.
				png_filename = paste0(prefix_str,fixed_vid_key, suffix_str)  
				png_file = list.files(pattern=png_filename, full.names=FALSE, recursive=TRUE)
				png_file_empty = vector.is.empty(png_file)
				
				if(png_file_empty==FALSE){
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
							 }
									
									
							  roi_polygons<-suppressMessages(left_join(ROIs_combined,df) %>% mutate(X_um = x_pixel/pixel_to_um, 
							                                          Y_um = y_pixel/pixel_to_um)
							  										)
							  roi_info<- suppressMessages(left_join(roi_info,df) )#%>% mutate(across(all_of(vars_to_plot), signif(.x,digits=3) ) )
								
							  upper_range_x = 256/pixel_to_um#ifelse(range(roi_polygons$x_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um
							  upper_range_y = 256/pixel_to_um#ifelse(range(roi_polygons$y_pixel)[2] > 256, range(roi_polygons$y_pixel)[2], 256)/pixel_to_um 
							  source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
						
							  ### merge stuff with real data here
							  na.value.forplot <- NA

							 				protocol_str = case_when(protocol_key == "singleAP" ~ "Single AP",
							 											protocol_key == "PP60" ~ "60 ms ISI",
							 											protocol_key == "PP75" ~ "75 ms ISI",
							 											protocol_key == "PP100" ~ "100 ms ISI",
							 											protocol_key == "PP150" ~ "150 ms ISI",
							 											protocol_key == "PP500" ~ "500 ms ISI")
							 				Ca_str = case_when(Ca_key == "0.5Ca" ~"0.5 mM Ca2+",
							 										Ca_key == "1Ca" ~ "1 mM Ca2+",
							 										Ca_key == "2Ca" ~ "2 mM Ca2+",
							 										Ca_key == "2Ca" ~ "4 mM Ca2+"
							 										)


							 				###### arrow math stuff ####
							 				ROI_ID_A<- roi_info %>% dplyr::filter(ROINumber == ROI_IDs[1]) %>% mutate(head_end_x = X_um,
							 																							head_end_y = Y_um - 0.2,
							 																							tail_start_x = X_um,
							 																							tail_start_y = Y_um -2)
							 				ROI_ID_B<- roi_info %>% dplyr::filter(ROINumber == ROI_IDs[2]) %>% mutate(head_end_x = X_um,
							 																							head_end_y = Y_um - 0.2,
							 																							tail_start_x = X_um,
							 																							tail_start_y = Y_um -2)
							 				scalebar = data.frame(x_start = 15 ,
							 										x_end = 25,
							 										y_start = 25,
							 										y_end = 25)

							 				 # scalebar = data.frame(x_start = user_x,
							 				# 						x_end = user_x+10,
							 				# 						y_start = user_y,
							 				# 						y_end = user_y)





											

											 for(k in 1:length(vars_to_plot)) {
											 	var=vars_to_plot[k]
											 		roi_info_mutated <- roi_info #%>% mutate(PPR_adj = signif(PPR, digits=2) )


											 	#if(guide_bool == TRUE) {
														 	p2<- ggplot() +
															  annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
															  geom_polygon(data=roi_polygons, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",fill=var,colour=var),alpha=1,size=4)+
															  scale_fill_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=guide_limits,name="PPR", oob=scales::squish)+
															  scale_colour_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=guide_limits, oob=scales::squish,guide="none")+
															  #geom_segment(data=ROI_ID_A, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  							  #geom_segment(data=ROI_ID_B, aes(x = tail_start_x, y = tail_start_y, xend = head_end_x, yend = head_end_y),size=2,colour="white",arrow = arrow(length = unit(0.5, "cm")))+
								  							  geom_segment(data=scalebar, aes(x=x_start,y=y_start,xend=x_end,yend=y_end), size= 5, colour="white")+
								  

															  #scale_color_manual(values = 'black', labels = 'Missing value') +
											  				  #geom_text(data=roi_info_mutated, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",label="PPR_adj"),  colour="black", size=5, fontface='bold')+
															  
															  labs(x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), title = paste0(Ca_str,", ",protocol_str),tag=tag_var ) +  #," : ",fixed_vid_key)) +
															  theme_minimal()+
															  scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
															  scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+
															  my.theme_spatialROI+
															  guides(fill=guide_legend(title=guide_title,override.aes = list(size = 1.2)))+
															  theme(legend.justification = c(1, 1), legend.position = c(1, 1),
																			  		legend.box.background = element_rect(fill = "white", color = "black"),
																			  		axis.text = element_blank(),
																			  		axis.title=element_blank(),
																			  		legend.title = element_text(colour="black", size=22, family="sans")#,
																			  		#plot.title=element_blank()
                 																	)
															  

															  
															}
											print(paste0("We should be saving the following file: ",Ca_key,"_",protocol_key,"_",fixed_vid_key,"_",var,"_spatial_plot_.png") )
											ggsave(filename=paste0(Ca_key,"_",protocol_key,"_",fixed_vid_key,"_",var,"_spatial_plot_.png"),plot=p2, device="png",dpi=600, bg="white",units="in",width=16,height=16)
}






	
							  	# p1<- ggplot() +
								#   annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
								#   geom_polygon(data=roi_polygons, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber"),colour="black",fill="white",alpha=0.5,size=0.5)+
								#   #scale_fill_gradient(low="blue",high="orange",na.value=na.value.forplot,limits=c(0,3),alpha=1)+
								#   #scale_color_manual(values = 'black', labels = 'Missing value') +
				  				#   #geom_text(data=roi_info_mutated, aes_string(x=("X_um"), y=("Y_um"),group="ROINumber",label="PPR_adj"),  colour="black", size=4, fontface='bold')+
								  
								#   labs(x = expression("X ("*mu*"m)"), y = expression("Y ("*mu*"m)"), title = paste0(fixed_vid_key,": Spatial Plot with ROI overlay only")) +
								#   theme_minimal()+
								#   scale_y_reverse(expand=c(0,0), limits=c(upper_range_y,0), oob=scales::oob_keep)+
								#   scale_x_continuous(expand=c(0,0),limits=c(0,upper_range_x), oob=scales::oob_keep)+
								#   #guides(color = guide_legend(override.aes = list(fill = na.value.forplot)))+
								
								#   my.theme_spatialROI
								
								
								# 	ggsave(filename=paste0(fixed_vid_key,"_spatial_plot_onlyROIs.png"),plot=p1, device="png",dpi=600, bg="white",units="in",width=16,height=16)
								












												
															  
												# } else if(guide_bool == FALSE){
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
												

													# if(save_bool == TRUE) {
													# 	print("this imaging region should be saved")
													# 	#tmp_plot <<- p2
													# 	return(p2)
													# 	ggsave(filename=paste0(fixed_vid_key,"_spatial_plot_.png"),plot=p1, device="png",dpi=600, bg="white",units="in",width=16,height=16)
													# } else { print("something strange happened")
													
													#return(print("something strange happened"))
												#}

				#}

	#}


