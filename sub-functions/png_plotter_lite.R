



library(ggplot2)
library(grid)
library(dplyr)
library(stringr)
library(imager)
library(readr)
library(bioimagetools)
library(EBImage)
library(ggthemes)




png_plotter_lite<- function(png_dir = png_dir,png_name=png_name, tag_var = taggy){



				vector.is.empty <- function(x) return(length(x) ==0)  



				png_dir <- png_dir
				setwd(png_dir)
				setwd(png_dir)


				
				# Specify the file name you want to find
					

				# List all directories (including subdirectories) in the parent directory
				#all_dirs <- list.dirs(png_dir, full.names = TRUE)

				#Derive the .png filename to look for.
				png_filename = png_name  
				png_file = list.files(pattern=png_filename, full.names=FALSE, recursive=TRUE)
				
				png_file_empty = vector.is.empty(png_file)
				
				if(png_file_empty==FALSE){
								y <- readImage(png_file)
								p1<- ggplot() + annotation_raster(y, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
												labs(tag=tag_var)+
												theme_void()+
												theme(plot.tag = element_text(colour="black", size=36, family="sans",face="bold"))

						}

				return(p1)
}