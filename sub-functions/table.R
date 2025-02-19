subDir <- "tableOutput"   


source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 
source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/myTheme.R") 


library(ggforce)



getTable = function(df = peakStats, grouper1 = grouper1, grouper2=NULL, comparisons = NULL) {


				se <- function(x) sqrt(var(x) / length(x))
				optimal.scale <- function(w,h, wanted.w, wanted.h) max(c(w/wanted.w,  h/wanted.h))

				#getTitle<- grouper1
				if (is.null(grouper2)) {


					print("no grouper2 anywhere!")
					grouping_vars = c(grouper1)
					print(grouping_vars)
				table = suppressMessages(df %>% group_by_at(.vars = grouping_vars) %>% 
								mutate(plateSegment = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,neuronSegment),
										plateRegion = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,region,neuronSegment),
										plateRegionROIs = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,region,neuronSegment,ROINumber)) %>%
								summarise(unique_coverslips = length(unique(plateSegment)),
											unique_neuronSegments = length(unique(plateRegion)),
											unique_ROIs = length(unique(plateRegionROIs)), 
											group_observations = n(),
											ROI_region_ratio = round(unique_ROIs/unique_neuronSegments, digits=2),
											obs_ROI_ratio = round(group_observations/unique_ROIs,digits=2),
											mean_dFF = round(mean(amplitude),2),
											se_dFF = round(se(amplitude),2),
											mean_halfwidth_ms = round(mean(halfwidth_ms),2),
											se_halfwidth_ms = round(se(halfwidth_ms),2),
											mean_tau_decay_ms = round(mean(tau_decay_ms),2),
											se_tau_decay_ms = round(se(tau_decay_ms),2)
														) 
											)
				
				
				tg = gridExtra::tableGrob(table)
  				h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
  				w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
  				scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape
				


				tableOutput = arrangeGrob(top = paste0("Descriptive statistics for dataset :", grouping_vars),tableGrob(table))

				ggsave(file=paste0(grouping_vars,"_table.pdf"), device='pdf', 
						plot = tableOutput, dpi=300, units="mm", 
						width=279, height=210, scale=scale)
				#table = df %>% group_by_at(grouper1)


				} else {


				grouping_vars = c(grouper1,grouper2)
				groupTitle = paste(grouping_vars[1],grouping_vars[2], sep='-')
				print("Since there's two groupers, grouping_vars = ")
				print(grouping_vars)
				table = suppressMessages(df %>% group_by_at(.vars = grouping_vars) %>% 
								mutate(plateSegment = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,neuronSegment),
										plateRegion = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,region,neuronSegment),
										plateRegionROIs = paste0(testProtein,Thrombin,Thrombin_state,dish,plate,region,neuronSegment,ROINumber)) %>%
								summarise(unique_coverslips = length(unique(plateSegment)),
											unique_neuronSegments = length(unique(plateRegion)),
											unique_ROIs = length(unique(plateRegionROIs)), 
											group_observations = n(),
											ROI_region_ratio = round(unique_ROIs/unique_neuronSegments, digits=2),
											obs_ROI_ratio = round(group_observations/unique_ROIs,digits=2),
											mean_dFF = round(mean(amplitude),2),
											se_dFF = round(se(amplitude),2),
											mean_halfwidth_ms = round(mean(halfwidth_ms),2),
											se_halfwidth_ms = round(se(halfwidth_ms),2),
											mean_tau_decay_ms = round(mean(tau_decay_ms),2),
											se_tau_decay_ms = round(se(tau_decay_ms),2)
														) 
								)

				
				
				tg = gridExtra::tableGrob(table)
  				h = grid::convertHeight(sum(tg$heights), "mm", TRUE)
  				w = grid::convertWidth(sum(tg$widths), "mm", TRUE)
  				scale = optimal.scale(w,h, 279, 210) + 0.1 #A4 = 279 x 210 in landscape
				


				tableOutput = arrangeGrob(top = paste0("Descriptive statistics for dataset :", grouping_vars),tableGrob(table))

				ggsave(file=paste0(groupTitle,"_table.pdf"), device='pdf', 
						plot = tableOutput, dpi=300, units="mm", 
						width=279, height=210, scale=scale)

				
				}


		tableOutput


}