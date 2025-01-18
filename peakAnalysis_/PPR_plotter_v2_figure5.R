#PPR_plotter_v1.R
#the goal of this script is to output a set of plots for EVERY ROI in the dataset across each Ca2+, and each pairedPulse paradigm. 
#This script will take advantage of existing code from peakFinder.R to output plots in pre-defined configuration. 
#The configuration of the output should be as follows:
#### PPR plot full ####
# nrow = 



subDir <- "PPR_plotter_v2_figure5_vSFN"   


source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/setFigureDirectory.R") 
#source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/myTheme.R") 
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/sub-functions/def_stim_vlines.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/subtract_by_index.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_v2.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_v3_ROIoverlay.R")
source("Z:\\R Scripts/Sam current as of 7.28.22/Active scripts/pipeline v5 - generalized/peakAnalysis/png_plotter_v4.R")

library(ggforce)
library(scales)
library(ggeasy)

library(ggpubr)
library(gridExtra)
library(gridtext)
library(grid)

PPR_stats<- function(df,groupers,stimEpoch_groupers,plotBy,levels, color_override =NULL,keys,save_ROIs = ROIs_to_save, keep_PP = keep_PP, tag_vars= tag_var_list){
						  

#### establish the plot parameters for the faceted, tracked ROIs
                        vid_keys_to_save = unique(str_extract(ROIs_to_save, "dish\\d\\d-plate\\d\\d-region\\d\\d") ) 

                        interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)
                        rm_groupers = c("exposeNum", "ROI_key","stimEpoch")
                        groupers_avg = groupers[!(groupers %in% rm_groupers)]
                        groupers_PPR = groupers_avg[!(groupers_avg %in% c("stimEpoch"))]

                        lineplot_width = 16
                        #lineplot_height = 24
                        scatterplot_width = 16
                        scatterplot_height = 16
                        color_switch = !is.null(color_override)
                        .keyvars = rlang::syms(keys)
                        rm_groupers = c("exposeNum")
                        CV_groupers = groupers[!(groupers %in% rm_groupers)] 

                        


                        df_subtracted<-subtract_by_index(df, groupers, groupers_avg)    
                        
                        PPR_calc<- df_subtracted%>% group_by_at(groupers_PPR) %>% dplyr::filter(#!is.na(stimEpoch), 
                                                                                                protocol != "singleAP", avg_normTime<0.85,avg_normTime>-0.85 ) %>%
                                                             summarise(interStim = (as.numeric( gsub( "PP", "", as.matrix( unique(protocol) ) ) ) / 1000),
                                                                        interStim_s = ifelse(!is.na(interStim), interStim, 0.5),
                                                                        time_of_stim1 = avg_normTime[which(reindex==0)],
                                                                        time_of_stim2 = time_of_stim1+interStim_s,
                                                                        cutoff_time2 = time_of_stim2+interStim_s,
                                                                        peak1 = max(avg_dFF[which( (avg_normTime > (time_of_stim1) ) & (avg_normTime < (time_of_stim2)) ) ], na.rm=TRUE),
                                                                        peak1_fix = ifelse(peak1<0.10, NA, peak1),
                                                                        peak2 = max(dFF_subtracted[which( (avg_normTime > (time_of_stim2) ) & (avg_normTime < (cutoff_time2) ) ) ], na.rm=TRUE), 
                                                                        peak2_fix = ifelse(peak2<0.10, NA, peak2), 
                                                                        PPR =ifelse(is.numeric(peak2_fix/peak1_fix), peak2_fix/peak1_fix, NA) 
                                                                        ) %>%
                                                                        mutate(pairedPulse_ISI_ms = as.numeric( gsub( "PP", "", as.matrix( protocol ) ) ) ) 
                        PPR_removedup = unique(PPR_calc)
                                                                        
                                                          


         PPRplot_data<- PPR_removedup %>%  dplyr::filter(Ca != "0pt25Ca")


        Ca_labels = c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'))

remove_labels <- function(original_func, remove_list = list()) {
  function(x) {
    original_result <- original_func(x)
    replace(original_result, original_result %in% remove_list, '')
  }
}


Ca.labs <- c("0.5 mM Ca2+", "1 mM Ca2+", "2 mM Ca2+")
names(Ca.labs) <- c("0pt5Ca", "1Ca","2Ca")
protocol.labs <- c("75 ms","100 ms", "500 ms")
names(protocol.labs) <- c("PP75", "PP100","PP500")



####### Title of Figure 5 (Figure 4?): Single synapses engage in diverse glutamate release behavior.    

######## A.  ###########
######## Protocol approach.                                               ###########


######## B. Plotting 3 trials for each protocol of an individual ROI, along with the df_subtracted version of the averaged trial. #########

######## C. Data with subtraction of the same ROI as in B. Play with different visualizations here. ########

######## D. Violin plots of Peak1 and Peak2 amplitudes ######
######## E. Violin plots of Peak2/Peak1 ######
########






##### AVG PPR WAVEFORMS ######
tmp_df <- df_subtracted %>% dplyr::filter(protocol %in% keep_PP)#, Ca %in% keep_Ca) 
tmp_df$protocol <- factor(tmp_df$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))

                    avg_PP_tracePlot<-ggplot(tmp_df, aes(x=avg_normTime, y = avg_dFF,  colour=Ca))+
                            geom_path(colour="grey61",size=0.5,alpha=0.2)+
                            stat_summary_bin(aes(group=Ca), geom="smooth", fun=mean, binwidth=interFrame.var*2,size=3)+
                            #stat_summary_bin(aes(group=Ca), geom="smooth", fun = mean, binwidth=interFrame.var,size=2)+
                            
                            #geom_vline(xintercept = 0, lty='dashed',colour='black',size=2)+
                         labs( x="Time (s)",
                                  y=expression(Delta*"F/F"),
                                  tag=tag_vars[1])+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=Ca_labels, values=color_override)}+
                            scale_x_continuous(breaks=c(-0.5,0,0.5),limits=c(-0.85,0.85))+

                            scale_y_continuous(breaks=c(0,0.5,1.0,1.5,2.0,2.5,3.0))+
                            theme_tufte()+
                            guides(colour="none",fill="none")+
                            coord_cartesian(ylim=c(-0.1,3.0),xlim=c(-1.1,1.1))+
                            my.theme+
                            facet_grid(protocol~Ca, labeller = labeller(protocol = protocol.labs, Ca = Ca.labs))

                            avg_peaks  <<- avg_PP_tracePlot
                                    rm(avg_PP_tracePlot)
                                
                
                            #{if(x_axis_switch==TRUE)easy_remove_x_axis()}


##### FULL PPRs #######

          PPRplot<-ggplot(PPRplot_data, aes(x=pairedPulse_ISI_ms, y=PPR,colour=Ca))+
                                             geom_hline(yintercept=1.0, lty="dashed", linewidth=2,colour="black",alpha=0.5)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="line", fun.y = mean,size=2,alpha=1)+
                                             #geom_smooth(aes(group=Ca,colour=Ca),formula=y~x, method=y~exp(x),alpha=1,se=FALSE,size=3)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="errorbar", fun.data=mean_se,width=20,size=2, alpha=1)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=6,alpha=1)+
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    title="",
                                                    tag = tag_vars[2]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=Ca_labels, values=color_override)}+
                                            coord_cartesian(ylim=c(0.5,1.5),xlim=c(50,550))+
                                            scale_y_continuous(breaks=c(0.5,1,1.5))+
                                            scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)), guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position = "none")
                                            
                                    avg_PPR  <<- PPRplot
                                    rm(PPRplot)
                                
                
                       

##### GET PPR VIOLINS

                            PPRplot_data$protocol <- factor(PPRplot_data$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))
    
                                   PPR_violin<-ggplot(PPRplot_data, aes(x=protocol,y=PPR,colour=Ca, fill=Ca))+
                                             geom_hline(yintercept=1.0, lty="dashed", linewidth=2,colour="black",alpha=0.5)+
                                            geom_violin(fill=NA, size=1.5,draw_quantiles = c(0.25, 0.5, 0.75))+
                                            geom_sina(size=2,alpha=0.7)+
                                            
                                            
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    title="",
                                                    tag=tag_vars[3]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_fill_manual(labels=Ca_labels, values=color_override)}+
                                            {if(color_switch)scale_colour_manual(labels=Ca_labels, values=color_override,guide='none')}+
                                            coord_cartesian(ylim=c(0,3))+
                                            scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3))+
                                            scale_x_discrete(labels=c("PP60" = "60", "PP75" = "75", "PP100" = "100","PP150" = "150","PP500" = "500" ), guide=guide_axis(angle=45))+
                                            #scale_x_continuous(expand=c(0,0))+
                                            
                                            facet_grid(~Ca, labeller = labeller( Ca = Ca.labs))+
                                            theme_tufte()+
                                            my.theme+
                                            guides(colour = guide_legend(show = FALSE),
                                                    fill = guide_legend(show = FALSE) )+
                                            theme(strip.text.y = element_blank(),
                                                    axis.text.x = element_text(size=36),
                                                    legend.position="none",
                                                    strip.text.x=element_blank())
                            
                                    violin_PPR  <<- PPR_violin

                                    rm(PPR_violin)




##### GET ALL THE SPATIAL DATA ##### 



###### for loop for plotting spatial data
        png_dir = "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/imageFiles/zstack/lowpwr_zstack/zstackOutputs_SQ"
        prefix_str = "^MAX_GluSnFR3_SyPhy_"
        suffix_str = "_0pt5Ca_zstack_lowpwr__FLATTENED.png"
        all_lut_str = "_0pt5Ca_zstack_lowpwr__FLATTENED_COLOR.png"
        dir_regex = "GluSnFR3(.*?)region\\d\\d_"

        subset_df<- PPRplot_data %>% dplyr::filter(vid_key %in% vid_keys_to_save,protocol %in% keep_PP)
        #print(unique(subset_df$protocol))
        vid_key_vec <- unique(subset_df$vid_key)
        Ca_vec <- unique(subset_df$Ca)
        protocol_vec <- unique(subset_df$protocol[which(subset_df$protocol != "singleAP")])
        #print(length(protocol_vec)) 
        vars_to_plot<- c("PPR")
        save_bool=TRUE
        levs <- c("singleAP","PP60","PP75","PP100","PP150","PP500")
        ordered_protocols <- protocol_vec[order(match(protocol_vec, levs))]
        print(ordered_protocols)

        tag_var = c("F","G","H",
                    "I","J","K",
                    "L","M","N")
        tag_var_subset = tag_vars[1:9]

        print("getting ROInames") 
        ROInames <- str_extract(save_ROIs,"ROI\\d\\d\\d\\d")   
        print(ROInames)
        for (k in 1:length(vid_key_vec)){
            
                vid_df<- PPRplot_data %>% dplyr::filter(vid_key == vid_key_vec[k])
                max_gradient_val = max(vid_df$PPR)

                
                    #if(unique(vid_df$vid_key) %in% vid_keys_to_save) { 
                    #    save_bool = TRUE
                    #    print(paste0( unique(vid_df$vid_key)," is a region we want to save out. Updated the save_bool value to TRUE" ) )
                    #} else {
                    #    save_bool = FALSE
                    #    #print("This is not a region we want to save out. Updated the save_bool value to FALSE")
                    #}

                    if(!is.null(vid_df) & save_bool==TRUE){
                    
                        png_plotter_v4(df = vid_df, png_dir = png_dir,prefix_str=prefix_str, suffix_str = suffix_str,all_lut_str = all_lut_str, dir_regex = dir_regex, vars_to_plot=vars_to_plot,save_bool=save_bool,tag_var = tag_vars[4])

                        png_plotter_ROIoverlay(df=vid_df,png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,vars_to_plot=vars_to_plot,save_bool=save_bool,tag_var = tag_vars[5],ROI_IDs = ROInames)
                                
                    }     

                count=0
                for(i in 1:length(Ca_vec)) {
                                    
                                    tmp_plot = NULL
                                    Ca_df<- vid_df %>% dplyr::filter(Ca == Ca_vec[i])
                                    check_data<- Ca_df %>% select(Ca, protocol, PPR)
                                    #print(check_data)
                                            
                                                for (j in 1:length(ordered_protocols)) {
                                                    protocol_df <- Ca_df %>% dplyr::filter(protocol == ordered_protocols[j])
                                                    if(length(protocol_df$protocol) > 0){

                                                    if(save_bool==TRUE) { 
                                                        #save_bool = TRUE
                                                        #print("This is not a protocol we want to save out. Updated the save_bool value to FALSE")
                                                            count = count+1
                                                            print(paste0(count, " was value of iteration counter at: Ca = ", unique(Ca_df$Ca), " and protocol = ", unique(protocol_df$protocol) ) )
                                                            if(count == 3) {
                                                               guide_bool = TRUE 
                                                            } else { guide_bool = FALSE}
                                                            tmp_tag_var = tag_var_subset[count]
                                                            tmp_plot = png_plotter(df=protocol_df,png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,vars_to_plot=vars_to_plot, save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag_var,ROI_IDs = ROInames)
                                                            #spat_facet_list
                                                            
                                                            nam <- paste("spat_facet", count, sep = "")
                                                            assign(nam, tmp_plot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                                                            }
                                                        

                                                        }

                                                    }
                                    
                                        }

                                }

        

                        
      



###### GENERATE SPECIFIC ROI DATA ####### 


  # Establish a tmp_data vector that can be used plot only the ROIs from a specific region for each iteration of the outer loop. 


PPR_data = suppressMessages(left_join(df_subtracted, PPR_removedup) %>% dplyr::filter(Ca != "0pt25Ca"))
tmp_PPR_data<- PPR_data %>% dplyr::filter(trackROI_key %in% save_ROIs)                           
tracked_ROIs <- unique(tmp_PPR_data$trackROI_key) 
Ca_labels = c(#expression("0.25 mM "*Ca^'2+'),
                    expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'), expression("4 mM "*Ca^'2+') )

for (i in 1:length(tracked_ROIs)) {
     
   print(paste0("loop iter on count: ", i))
   #print("Generating data visualizations for trackROI_key:", ROI)

   tag_var_subset = tag_vars[c(11,13)]
   tmp_data <- PPR_data %>%
         dplyr::filter(trackROI_key == tracked_ROIs[i])

    protocols = unique(tmp_data$protocol)

          tmp_PPRplot<-ggplot(tmp_data, aes(x=pairedPulse_ISI_ms, y=PPR,colour=Ca))+
                                            geom_hline(yintercept=1.0, lty="dashed", linewidth=2,colour="black",alpha=0.5)+
                                            
                                            geom_point(size = 3, alpha=0.4)+
                                            stat_summary(aes(group=Ca,colour=Ca),geom="line", fun.y=mean,alpha=1,size=3,na.rm=TRUE)+
                                            #stat_summary(aes(group=Ca, colour=Ca), geom="errorbar", fun.data=mean_se,width=20,size=2, alpha=1,na.rm=TRUE)+
                                            stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=6,alpha=1,na.rm=TRUE)+
                                             
                                            #stat_summary(aes(group=Ca, colour=Ca), geom="line", fun.y = mean,size=2,alpha=1)+
                                            #geom_point(size=5)+
                                            #geom_line(size=2.5)+
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    title=unique(tmp_data$ROINumber),
                                                    tag = tag_var_subset[i]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(labels=Ca_labels, values=color_override)}+
                                            coord_cartesian(ylim=c(0.5,2.5),xlim=c(50,550))+ 
                                            scale_y_continuous(breaks=c(0.5,1,1.5,2.0,2.5))+
                                            scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)))+
                            
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position="none")
                                            #theme(legend.position = c(.95, .95),
                                            #        legend.justification = c("right", "top"),
                                            #        legend.box.just = "right",
                                            #        legend.margin = margin(6, 6, 6, 6))
                                            
                       
                             nam <- paste("ROI_PPR_facet", i, sep = "")
                             assign(nam, tmp_PPRplot,envir = .GlobalEnv)

                         #  ggsave(filename=paste0("PPRplot_",ROI,"_.png"),plot=PPRplot, device="png",dpi=600, units="in",width=scatterplot_width,height=scatterplot_width)



    # Generate a faceted plot of each protocol
    
    levs <- c("singleAP","PP60","PP75","PP100","PP150","PP500")
    ordered_protocols <- protocols[order(match(protocols, levs))]
    #sort(ordered_protocols)
    #xmin = min(PPR_data$normTime,na.rm=TRUE)
    #xmax = max(PPR_data$normTime, na.rm=TRUE)
    #ymax = max(PPR_data$dFF)
    #ymin = min(PPR_data$dFF)
                    

    tmp = list()
    #for (j in 1:length( ordered_protocols) ) {

         library(tidyverse) 
         lineplot_width = 16
         lineplot_height = 20
         #tag_var = c("O","Q")
         tag_var_subset = tag_vars[c(10,12)]
         
       

        j_tmp_data <- tmp_data %>% dplyr::filter(protocol %in% keep_PP) %>% mutate(interStim = as.numeric( gsub( "PP", "", as.matrix( protocol ) ) ) / 1000 )#convert PP\\d... to \\d... and then to seconds)
        j_tmp_data$protocol <- factor(j_tmp_data$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))

                    x_int_data<- j_tmp_data %>% ungroup() %>% group_by(protocol) %>% summarise(x_int = mean(interStim))
                 
                    
                        tmp_tracePlot<-ggplot(j_tmp_data, aes(x=avg_normTime, y = avg_dFF,  colour=Ca)) +
#                             {if(vline_switch)geom_vline(xintercept = stimParadigm_times, colour = "red", alpha = 0.4, linetype = "solid", size = .9)} +
                            geom_path(alpha=1,size=3)+
                            #ggplot(j_tmp_data, aes(x=normTime, y = dFF)) +
                            #geom_vline(data=x_int_data,aes(xintercept = c(0,x_int),group=protocol), colour = "red", alpha = 0.4, linetype = "solid", size = .9) +
                            #geom_path(colour="grey21",alpha=0.5,size=2)+
                            #geom_line(aes(x=normTime, y=avg_dFF, colour=Ca),size=2)+
                            #stat_summary_bin(aes(group=Ca), geom="smooth", fun = mean, binwidth=interFrame.var,size=2)+
                            
                            labs( x="Time (s)",
                                  y=expression(Delta*"F/F"),
                                  title=unique(j_tmp_data$ROINumber),
                                  tag=tag_var_subset[i])+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=Ca_labels, values=color_override)}+
                            scale_x_continuous(breaks=c(-0.5,0,0.5),limits=c(-0.85,0.85))+

                            scale_y_continuous(breaks=c(0,1.0,2.0))+
                            theme_tufte()+
                            guides(colour="none",fill="none")+
                            coord_cartesian(ylim=c(-0.1,2.0),xlim=c(-1.1,1.1))+
                            #{if(axis_switch==TRUE)easy_remove_y_axis()}+

                            my.theme+
                            facet_grid(protocol~Ca, labeller = labeller(protocol = protocol.labs, Ca = Ca.labs))
                      
                            nam <- paste("ROI_PPR_facet_traces_",i, sep = "")
                            assign(nam, tmp_tracePlot,envir = .GlobalEnv)
              #          }



                }





    PPR_removedup

}						
						