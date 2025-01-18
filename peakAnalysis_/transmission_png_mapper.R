			##### pseudocode for Extended Figure 6


			###first, calculate the variables needed for viz
			### z-score ratio is prolly the ticket
			### then, generate arrows that can be viz
			### to generate map pre-treatment and map treatment

			### dplyr::filter(chemical_condition == "ctrl")
			### dplyr::filter(chemical_condition != "ctrl")


			#### function(ROI_to_save, tag_var, plot_var, guide_limits, guide_title, scale_bar_bool, scale_bar_x,scale_bar_y)





# In principle, transmission_png_map should generate a single png output each time it is called.
### transmission_png_map is a standalone container function which:
#1. accepts a dataframe containing summary data, 
#2. a string specifying an ROI, 
#3. a string specifying a folder name where .png files are stored,
#4. var_to_plot
#5. guide_title
#6. guide_lim
#7. tag_var








subDir <- "transmission_png_map"   


source(paste0(path,"sub-functions/setFigureDirectory.R")) 
source(paste0(path,"sub-functions/myTheme.R") )
source(paste0(path,"peakAnalysis/png_plotter_v5_ROIoverlay_syntransmission.R"))
source(paste0(path,"peakAnalysis/plot_dualRecords_fxn.R"))
source(paste0(path,"peakAnalysis/save_plot_as_jpeg_and_emf.R"))
source(paste0(path,"peakAnalysis/png_plotter_v2_synTrans.R"))



library(ggforce)
library(scales)
library(ggthemes)
library(ggpmisc)
library(ggpubr)


transmission_png_map<- function(summary_df, ROI_to_save, png_folder,var_to_plot,guide_title,tag_var){

                    vid_keys_to_save = unique(str_extract(save_ROIs, "dish(.*?)region\\d\\d") )
                    print(paste0("Vid_keys_to_save = ", vid_keys_to_save ))
                    if (!is.null(interSpike_thresh)) { interSpike_thresh = interSpike_thresh
                                                        } else {
                                                            interSpike_thresh = 0.5
                                                            } ## set threshold for accepting interSpike interval as transmission
                    secondAxis_switch = !is.null(secondAxis)
                    color_switch = !is.null(color_override)
                    comparison_switch = !is.null(comparisons)
                    ## generate the keys for tracking ROIs
                    barplot_width = 16
                    lineplot_width = 16
                    plot_height = 16
                    plot_dim=16




###### GET SPATIAL DATA ##### ROI OVERLAY + ARROW

# ###### for loop for plotting spatial data
        png_dirs = c(#"Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_3",
                        "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_6",
                        "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_11")
        prefix_str = "AVG_GluSnFR3_JF646_4Ca_cis_"
        suffix_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        all_lut_str = "_ctrl_repl01__FLATTENED_COLOR.png"
        dir_regex = "registered_(.*?)_"
        save_bool=TRUE
        #print(paste0("Available vid_keys to save out are: ", unique(traces_df$vid_key) ) ) 
        #subset_df<- traces_df_fixed %>% dplyr::filter(vid_key %in% vid_keys_to_save, trackROI_key %in% clean_ROI_list)
        #print(paste0("Column names are: ", colnames(subset_df)) ) 
        #print(paste0("Vid_key is: ", unique(subset_df$vid_key) ) )
        #print(unique(subset_df$protocol))
        #vid_key_vec <- unique(subset_df$vid_key)
        
        print("Taking a crack at generating png plots")
        subset_df<- transmission_rateCheck %>% dplyr::filter(vid_key %in% vid_keys_to_save,chemical_condition == "Pre-Treatment")

        #ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d") 
        vid_key_vec <- unique(subset_df$vid_key)
        vars_to_plot<- c('med_Glu_amplitude', "med_JF_amplitude",'JF_freq','transmission_rate')#c("releaseProbability_perROI")
        save_bool=TRUE
        guide_bool = TRUE

        guide_limits = list(c(0.5,2),c(0.1,0.8),c(0,0.25),c(0,1))
        guide_titles = c(expression(iGlu[Delta*F/F]),expression(JF[646][Delta*F/F]),expression(JF[646]~(Hz)),expression(P[transmission]))# # expression(CV[tau*"decay"])#
        
        #expression("Avg."~Evoked[Delta*"F/F"] / Spont[Delta*"F/F"])
        tmp_tag_var_override = c("F","G", "H")


        ###scale bar params! set the left-corner 
        user_x = 51.2-10.6
        user_y = 50.2
        tag_vars = c("F","H"    )



        tag_vars = c("H","F")
        print("getting ROInames") 
        ROInames <- str_extract(save_ROIs,"ROI\\d\\d\\d\\d")   
        print(ROInames)
        for(j in 1:length(png_dirs)){
            current_png_dir = png_dirs[j]
            print(paste0("current png directory is: ",current_png_dir) )
            print(paste0("Vid_key is : ",vid_key_vec[j]) )
               
            tmp_subset_transmission <- subset_df %>% dplyr::filter(vid_key == vid_key_vec[j]) %>% mutate(fix_Ca = "4Ca")
            print(head(tmp_subset_transmission))

            for(k in 1:length(vars_to_plot)){            
                                guide_bool = TRUE
                                var_to_plot = vars_to_plot[k]
                                tmp_guide_title = guide_titles[k]
                                tmp_limits = guide_limits[[k]]
                                tmp_tag = ""#tmp_tag_var_override[j]
                                if(k == 1){tmp_scalebar = TRUE} else {tmp_scalebar = FALSE}
                                # tmp_plot = png_plotter(df=tmp_subset_transmission,
                                #                                     png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                #                                     vars_to_plot=var_to_plot, 
                                #                                     save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
                                #                                     ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar,add_arrows=FALSE,
                                #                                     scalebar_length=10,binfactor=2)

            }
        }

        for (k in 1:total_ROIs){
                subset_traces = traces_df_fixed %>% dplyr::filter(trackROI_key %in% clean_transmission_ROIs_list[k])
                current_ROIs = unique(subset_traces$ROINumber)            
                
               current_png_dir = png_dirs[k]
               current_tag = tag_vars[k]
               save_bool=TRUE    
                    
                       png_plot<-png_plotter_ROIoverlay(df=subset_traces,png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,save_bool=save_bool,tag_var = current_tag,ROI_IDs = current_ROIs)
                        

                        if(k == 1){
#### save to global env
                            regionH  <<- png_plot


                            }
                        if(i == 2){
#### save to global env
                           regionF  <<- png_plot

                            }         
                         
                    

        }