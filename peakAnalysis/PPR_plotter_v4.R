#PPR_plotter_v1.R
#the goal of this script is to output a set of plots for EVERY ROI in the dataset across each Ca2+, and each pairedPulse paradigm. 
#This script will take advantage of existing code from peakFinder.R to output plots in pre-defined configuration. 
#The configuration of the output should be as follows:
#### PPR plot full ####
# nrow = 



subDir <- "PPR_plotter_outputs"   


source(paste0(path, "sub-functions/setFigureDirectory.R")) 
#source(paste0(path, "sub-functions/myTheme.R") 
source(paste0(path, "sub-functions/myTheme.R") )
source(paste0(path, "sub-functions/def_stim_vlines.R"))
source(paste0(path, "peakAnalysis/subtract_by_index.R"))
source(paste0(path, "peakAnalysis/png_plotter_v2_basalRP.R"))
source(paste0(path, "peakAnalysis/save_plot_as_jpeg_and_emf.R"))

library(ggforce)
library(scales)
library(ggeasy)

library(ggpubr)
library(gridExtra)
library(gridtext)
library(grid)

PPR_stats<- function(df,groupers,stimEpoch_groupers,plotBy,levels, color_override =NULL,keys,save_ROIs = ROIs_to_save, keep_PP = keep_PP, tag_vars= tag_var_list, png_dir = png_dir){
						  

#### establish the plot parameters for the faceted, tracked ROIs
                        vid_keys_to_save = unique(str_extract(ROIs_to_save, "dish\\d\\d-plate\\d\\d-region\\d\\d") ) 

                        interFrame.var = as.vector(df %>% ungroup() %>% summarise(interFrame=mean(interFrame,na.rm=TRUE))) %>% unlist(., use.names=FALSE)
                        print(paste0("Interframe.var = ", interFrame.var))


                        check<- df %>% dplyr::filter(normTime >= 0) %>% group_by(ROI_key) %>% summarise(max_frameTime = max(interFrame))
                        ROI_keys_to_omit <- check %>% dplyr::filter(max_frameTime > 0.18)
                        ROI_keys_to_omit <- unique(ROI_keys_to_omit$ROI_key)



                        df_fixed_width <- df %>% dplyr::filter(normTime > -0.5,normTime < 1.00) %>% dplyr::filter(!ROI_key %in% ROI_keys_to_omit ) %>% dplyr::filter(Ca != "0pt5Ca")
                        # maxTime = max(df_fixed_width$normTime)
                        # minTime = min(df_fixed_width$normTime)
                        # print(paste0("Max time is = ", maxTime))
                        # print(paste0("min time is = ", minTime))
                        # bins_for_ntile = ceiling( (maxTime - minTime) / interFrame.var )
                        # print(paste0("Number of bins for ntile is = ", bins_for_ntile))
                        
                        rm_groupers = c("exposeNum", "ROI_key","stimEpoch")
                        groupers_avg = groupers[!(groupers %in% rm_groupers)]
                        groupers_PPR = groupers_avg[!(groupers_avg %in% c("stimEpoch"))]

                        plot_height = 16 #16
                        #lineplot_height = 24
                        color_switch = !is.null(color_override)
                        .keyvars = rlang::syms(keys)
                        rm_groupers = c("exposeNum")
                        CV_groupers = groupers[!(groupers %in% rm_groupers)] 

                        

                        df_subtracted<-subtract_by_index(df=df_fixed_width,groupers=groupers,groupers_avg=groupers_avg,ntile_bins = bins_for_ntile)    
                        
                        # PPR_calc<- df_subtracted%>% group_by_at(groupers_PPR) %>% dplyr::filter(#!is.na(stimEpoch), 
                        #                                                                         protocol != "singleAP", avg_normTime<0.85,avg_normTime>-0.85 ) %>%
                        #                                      summarise(interStim = (as.numeric( gsub( "PP", "", as.matrix( unique(protocol) ) ) ) / 1000),
                        #                                                 interStim_s = ifelse(!is.na(interStim), interStim, 0.5),
                        #                                                 time_of_stim1 = avg_normTime[which(reindex==0)],
                        #                                                 time_of_stim2 = time_of_stim1+interStim_s,
                        #                                                 cutoff_time2 = time_of_stim2+interStim_s,
                        #                                                 peak1 = max(avg_dFF[which( (avg_normTime > (time_of_stim1) ) & (avg_normTime < (time_of_stim2)) ) ], na.rm=TRUE),
                        #                                                 peak1_fix = ifelse(peak1<0.10, NA, peak1),
                        #                                                 peak2 = max(dFF_subtracted[which( (avg_normTime > (time_of_stim2) ) & (avg_normTime < (cutoff_time2) ) ) ], na.rm=TRUE), 
                        #                                                 peak2_fix = ifelse(peak2<0.10, NA, peak2), 
                        #                                                 PPR =ifelse(is.numeric(peak2_fix/peak1_fix), peak2_fix/peak1_fix, NA) 
                        #                                                 ) %>%
                        #                                                 mutate(pairedPulse_ISI_ms = as.numeric( gsub( "PP", "", as.matrix( protocol ) ) ) ) 


                         PPR_calc<- df_subtracted%>% group_by_at(groupers) %>% dplyr::filter(#!is.na(stimEpoch), 
                                                                                                protocol != "singleAP", normTime<0.85,normTime>-0.85 ) %>%
                                                             summarise(interStim = (as.numeric( gsub( "PP", "", as.matrix( unique(protocol) ) ) ) / 1000),
                                                                        interStim_s = ifelse(!is.na(interStim), interStim, 0.5),
                                                                        time_of_stim1 = normTime[which(reindex==0)],
                                                                        time_of_stim2 = time_of_stim1+interStim_s,
                                                                        cutoff_time2 = time_of_stim2+interStim_s,
                                                                        peak1 = max(dFF[which( (normTime > (time_of_stim1) ) & (normTime < (time_of_stim2)) ) ], na.rm=TRUE),
                                                                        time_of_peak1 = normTime[which(dFF == peak1)],
                                                                        peak1_fix = ifelse(peak1<0.35, NA, peak1),
                                                                        peak2 = max(dFF_subtracted[which( (normTime > (time_of_stim2) ) & (normTime < (cutoff_time2) ) ) ], na.rm=TRUE), 
                                                                        time_of_peak2 = normTime[which(dFF_subtracted == peak2)],
                                                                        peak2_fix = ifelse(peak2<0.35, NA, peak2), 
                                                                        PPR =ifelse(is.numeric(peak2_fix/peak1_fix), peak2_fix/peak1_fix, NA) 
                                                                        ) %>%
                                                                        mutate(pairedPulse_ISI_ms = as.numeric( gsub( "PP", "", as.matrix( protocol ) ) ) ) 
                        PPR_removedup = unique(PPR_calc)


                        print(paste0("We've quantified this many synapses:", length(unique(PPR_removedup$trackROI_key))))

                        PPR_summary = PPR_removedup %>% group_by(Ca, Ca_mM, protocol,vid_key,trackROI_key,ROINumber,pairedPulse_ISI_ms) %>% summarise(n_obs = n(),
                                                                                                                                 Ca_expr = case_when(Ca_mM == 0.5 ~"0pt5Ca",
                                                                                                                                                    Ca_mM == 1 ~"1Ca",
                                                                                                                                                    Ca_mM == 2 ~"2Ca"),
                                                                                                                                sum_PPR_detect_na = sum(is.na(PPR)),
                                                                                                                                mean_PPR_perROI = mean(PPR,na.rm=TRUE),
                                                                                                                                PPR_category = case_when(#sum_PPR_detect_na > 0 & mean_PPR_perROI == NaN ~ "Secret Non-responder",
                                                                                                                                                         is.na(mean_PPR_perROI) ~ "Non-responder",
                                                                                                                                                                #mean_PPR_perROI < 0.5 ~ "Strong depressing",
                                                                                                                                                                mean_PPR_perROI >= 1.6 ~ "Strong facilitating",
                                                                                                                                                                mean_PPR_perROI >= 1.1 & mean_PPR_perROI < 1.6 ~ "Facilitating",
                                                                                                                                                                mean_PPR_perROI >= 0.9 & mean_PPR_perROI < 1.1 ~ "Neutral",
                                                                                                                                                                #mean_PPR_perROI >= 1.0 & mean_PPR_perROI < 1.5 ~ "Weak facilitating",
                                                                                                                                                                mean_PPR_perROI < 0.9 ~ "Depressing")
                                                                                                                                )
                                                                        
                        plot_avgs = PPR_summary %>% group_by(Ca, Ca_mM, protocol, pairedPulse_ISI_ms) %>% summarise(n_all = sum(n_obs),
                                                                                                                        mean_PPR = mean(mean_PPR_perROI,na.rm=TRUE),
                                                                                                                        sd_PPR = sd(mean_PPR_perROI,na.rm=TRUE),
                                                                                                                        se_PPR = sd_PPR/sqrt(n_all)
                                                                                                                        )                           


         PPRplot_data<- PPR_removedup# %>%  dplyr::filter(Ca != "0pt25Ca")




Ca.labs <- c(#expression("0.5 mM "*Ca^'2+'), 
            expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'))  

#names(Ca.labs) <- c("0pt5Ca", "1Ca","2Ca")

protocol.labs <- c("Test~Pulse", "60~ms", "75~ms","100~ms","150~ms", "500~ms")

remove_labels <- function(original_func, remove_list = list()) {
  function(x) {
    original_result <- original_func(x)
    replace(original_result, original_result %in% remove_list, '')
  }
}



cc <- scales::seq_gradient_pal("blue", "orange", "Lab")(seq(0,1,length.out=4))
parse.labels <- function(x) parse(text = x)
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_blank()#element_text(size=14, face="bold")
  )

cc<- c("grey1",cc)
pie_chart_data<- PPR_summary %>% ungroup() %>% group_by(Ca, Ca_mM,protocol,PPR_category) %>% summarise(counts = length(unique(trackROI_key) ) )
pie_chart_total<- pie_chart_data %>% group_by(Ca,Ca_mM,protocol) %>% summarise(total_counts = sum(counts))

pie_chart_final<- left_join(pie_chart_data, pie_chart_total) %>% mutate(proportion = counts/total_counts,
                                                                        asPercent = proportion*100)

PPR_bias_for_spatial_plot<- PPR_summary %>% 
                        ungroup() %>% group_by(Ca, Ca_mM,trackROI_key) %>%
                        dplyr::filter(Ca_mM == 1) %>%  
                        dplyr::filter(!PPR_category %in% c("Neutral", "Depressing",  "Non-responder")) %>%#c("Neutral", "Depressing",  "Non-responder")) %>%
                        mutate(max_PPR = max(mean_PPR_perROI,na.rm=TRUE),
                                PPR_bias = case_when(mean_PPR_perROI == max_PPR ~ protocol,
                                                     mean_PPR_perROI != max_PPR ~ NA),
                                PPR_bias_numeric = case_when(PPR_bias == "PP60" ~ "60 ms",
                                                             PPR_bias == "PP75" ~ "75 ms",
                                                             PPR_bias == "PP100" ~ "100 ms",
                                                             PPR_bias == "PP150" ~ "150 ms",
                                                             PPR_bias == "PP500" ~ "500 ms")
                                ) %>%
                        select(Ca,Ca_mM,vid_key,protocol,trackROI_key,ROINumber,mean_PPR_perROI,max_PPR,PPR_bias,PPR_bias_numeric)
PPR_bias_for_spatial_plot <- PPR_bias_for_spatial_plot %>% dplyr::filter(!is.na(PPR_bias))
           
print(PPR_bias_for_spatial_plot)

PPR_bias<- PPR_summary %>% 
                        ungroup() %>% group_by(Ca, Ca_mM,trackROI_key) %>%  
                        dplyr::filter(!PPR_category %in% c("Neutral","Depressing",  "Non-responder")) %>%
                        mutate(max_PPR = max(mean_PPR_perROI,na.rm=TRUE),
                                PPR_bias = case_when(mean_PPR_perROI == max_PPR ~ protocol,
                                                     mean_PPR_perROI != max_PPR ~ NA)
                                ) %>%
                        ungroup() %>%
                        dplyr::filter(!is.na(PPR_bias))%>%
                        group_by(Ca, Ca_mM,PPR_bias) %>%
                        summarise(counts =length(unique(trackROI_key)))
PPR_bias_total<- PPR_bias %>% group_by(Ca,Ca_mM) %>% summarise(total_counts = sum(counts))
PPR_bias_final <- left_join(PPR_bias, PPR_bias_total) %>% mutate(proportion = counts/total_counts,
                                                                        asPercent = proportion*100,
                                                                        pairedPulse_ISI_ms = as.numeric( gsub( "PP", "", as.matrix( PPR_bias ) ) ) 

                                                                        )



print(PPR_bias_final)


  barPlot <-ggplot(PPR_bias_final, aes(x=pairedPulse_ISI_ms,y=proportion, colour=Ca)) +
                                                
                                                geom_line(alpha=1,size=1.5,lty="solid")+#,colour='#22A884FF')+              
                            
                                                geom_point(shape=21,stroke=1.5,size=4,fill="white")+#,colour='#22A884FF')+
                                                
                                                labs( subtitle = "Facilitation Preference",
                                                      x=expression("ISI of Max PPR (ms)"),
                                                      y="Fraction of boutons",
                                                      tag = "D")+
                                                #scale_colour_manual())+
                                                {if(color_switch)scale_colour_manual(labels=Ca.labs, values=color_override)}+
                            
                                                scale_x_continuous(breaks=c(60,75,100,150,500),labels = remove_labels(scales::label_number_auto(), c(60,100)))+
                                                scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels = number_format(accuracy = 0.01))+
                                                coord_cartesian(ylim=c(0,0.6))+
                                                theme_tufte()+
                                                my.theme+
                                                theme(legend.position = "none")
                                                # theme(#legend.spacing.y = unit(2.0,'cm'),
                                                #        # legend.position = c(1, .5),
                                                #         #legend.justification = "right",
                                                #         #legend.box.just = "right",
                                                #         #legend.margin = margin(6, 6, 6, 6)
                                                #         )#+
                                                
                                            save_plot(filename=paste("synapseProportions_PPR_bias", sep='_'),plot=barPlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                            PPR_bias<<-barPlot




print(pie_chart_final)

pie_chart_final$PPR_category = factor(pie_chart_final$PPR_category, levels = c("Non-responder","Strong depressing","Depressing","Neutral", "Weak facilitating","Facilitating","Strong facilitating"))#("Inactive:~P[iGlu]~=~0","Rarely~active:~0.25~>~P[iGlu]~>~0", "Active:~0.75~>~P[iGlu]~>~0.25","Very~active:~P[iGlu]~>~0.75"))
pie_chart_final$Ca = factor(pie_chart_final$Ca, levels = c("1Ca","2Ca")) #"0pt5Ca",
pie_chart_final$protocol = factor(pie_chart_final$protocol, levels = c("PP60","PP75","PP100","PP150","PP500"))
levels(pie_chart_final$Ca) <- Ca.labs
levels(pie_chart_final$protocol) <- protocol.labs[-1]

#activity_labels = c(expression("Inactive: "*P[iGlu]~"="~0),expression("Rarely Active: "*0~"<"~P[iGlu]~"<"~0.25),expression("Somewhat Active: "*0.25~"\u2264"~P[iGlu]~"<"~0.50),expression("Active: "*0.50~"\u2264"~P[iGlu]~"<"~0.75),expression("Very Active: "*0.75~"\u2264"~P[iGlu]),expression("Multivesicular") )

#levels(pie_chart_final$activity) = activity_labels

pie<- ggplot(pie_chart_final, aes(x="", y=asPercent, fill=PPR_category))+
        geom_bar(width = 1, stat = "identity",colour="black",size=0.05)+
        coord_polar("y", start=0)+
        scale_fill_manual(values=cc)+
        #labs(tag="B")+
        guides(fill = guide_legend(reverse = TRUE,label.position="left"))+
        facet_grid(Ca~protocol, labeller=label_parsed)+#, ncol=5,nrow=3)+
        #my.theme+
        blank_theme +
        theme(strip.text = element_text(colour="black", size=scalefactor*22, family="sans"),
                plot.tag = element_text(colour="black", size=scalefactor*48, family="sans",face="bold"),
                legend.text =element_text(colour="black", size=scalefactor*22, family="sans"),
                legend.position="left",
                legend.justification = c("left","center"),
                axis.text = element_blank(),
                legend.title=element_blank())


            # pirid <- ggplot(pie, aes(x = "", y = Counts, fill = Puncta)) +
            #           geom_col(color = "black") +
            #           geom_text(aes(label = percent),
            #                     position = position_stack(vjust = 0.5)) +
            #           coord_polar(theta = "y") + theme(axis.ticks = element_blank(),
            #                 axis.title = element_blank(),
            #                 axis.text = element_text(size = 15), 
            #                 legend.position = "none", # Removes the legend
            #                 panel.background = element_rect(fill = "white")) + 
            #           scale_color_viridis(option = "D", discrete = TRUE) +
            #           #scale_fill_viridis(discrete = TRUE) +
            #           theme_void()#+
            #           #my.theme()
save_plot(filename=paste0("piechart_trial"),plot=pie, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)


####### Title of Figure 5 (Figure 4?): Single synapses engage in diverse glutamate release behavior.    

######## A.  ###########
######## Protocol approach.                                               ###########


######## B. Plotting 3 trials for each protocol of an individual ROI, along with the df_subtracted version of the averaged trial. #########

######## C. Data with subtraction of the same ROI as in B. Play with different visualizations here. ########

######## D. Violin plots of Peak1 and Peak2 amplitudes ######
######## E. Violin plots of Peak2/Peak1 ######
########



Ca.labs <- c(#expression("0.5 mM "*Ca^'2+'), 
              expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'))  

names(Ca.labs) <- c( "1Ca","2Ca") #"0pt5Ca",

protocol.labs <- c("Test~Pulse", "60~ms", "75~ms","100~ms","150~ms", "500~ms")

#names(protocol.labs) <- c("singleAP", "PP60", "PP75", "PP100","PP150","PP500")


#Ca_levels<- c(expression("0.5 mM "*Ca^'2+'), expression("1 mM "*Ca^'2+'), expression("2 mM "*Ca^'2+'))  
#protocol_levels =                         
                        
                        

##### AVG PPR WAVEFORMS ######
tmp_df <- df_subtracted %>% #dplyr::filter(protocol %in% keep_PP) %>% 
                                  mutate(stim2_vline = as.numeric(gsub("PP", "", protocol) )/1000,
                                          Ca_expr = case_when(Ca_mM == 0.5 ~"0pt5Ca",
                                                              Ca_mM == 1 ~"1Ca",
                                                              Ca_mM == 2 ~"2Ca"),
                                          protocol_expr = as.factor(as.character(protocol)) 
                                                                            )
tmp_df$Ca_expr <- factor(tmp_df$Ca_expr, levels = c("1Ca","2Ca")) #"0pt5Ca",
levels(tmp_df$Ca_expr) <- Ca.labs

tmp_df$protocol_expr <- factor(tmp_df$protocol_expr, levels = c("singleAP","PP60","PP75","PP100","PP150","PP500"))
levels(tmp_df$protocol_expr) <- protocol.labs



#levels(tmp_df$Ca) <- Ca_levels

#print(paste0(  "Max interframe for  1 mM Ca2+, 500 ms ISI, occurred at the following ISIs: ", unique(check$max_frameTime) ) )
#print(unique(check$vid_key))
vline_df <- tmp_df %>% select(Ca, protocol,Ca_expr,protocol_expr, stim2_vline)
vline_df <- unique(vline_df)

vline_df$protocol <-factor(vline_df$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))
tmp_df$protocol <- factor(tmp_df$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))
tmp_df$Ca<- factor(tmp_df$Ca, levels = c("1Ca","2Ca")) #"0pt5Ca",

# levels(tmp_df$Ca) = Ca.labs
# levels(tmp_df$protocol) = protocol.labs
# levels(vline_df$protocol) = protocol.labs

#                     avg_PP_tracePlot<-ggplot(tmp_df, aes(x=normTime, y = dFF,  colour=Ca))+
#                             #geom_hline(yintercept = c(1,2,3,4), colour='black',size=0.5,alpha=0.1,lty="solid")+
#                             #geom_path(colour="grey61",size=1,alpha=0.2)+
#                             #geom_smooth(method='loess')+
#                             stat_summary_bin(aes(group=Ca), geom="smooth", fun=mean, binwidth=interFrame.var*2,size=2.5,alpha=1)+
#                             #stat_summary_bin(aes(group=Ca), geom="smooth", fun = mean, binwidth=interFrame.var,size=2)+
                            
#                             geom_vline(xintercept = 0, lty='longdash',colour='black',size=1)+
#                             geom_vline(data=vline_df,aes(xintercept = stim2_vline, group=Ca), lty='longdash',colour='black',size=1)+
#                          labs( x="Normalized Time (s)",
#                                   y=expression(Delta*"F/F"),
#                                   tag="")+#tag_vars[1])+
#                             {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
#                             {if(color_switch)scale_colour_manual(labels=Ca.labs, values=color_override)}+
#                             scale_x_continuous(breaks=c(0,0.5),limits=c(-0.80,0.80))+

#                             scale_y_continuous(breaks=c(0,0.5,1.0,1.5,2.0,3.0,4.0))+
#                             theme_tufte()+
#                             guides(colour="none",fill="none")+
#                             coord_cartesian(ylim=c(0,1.5),xlim=c(-0.2,0.80))+
#                             my.theme+
#                             facet_grid(protocol_expr~Ca_expr, labeller = label_parsed )+ 
#                                                             #labeller(protocol = protocol.labs, Ca = as_labeller(Ca.labs,label_parsed))+
#                             theme(legend.position = "none",
#                                                   strip.text.y.right = element_text(angle = 0))

#                             #avg_peaks  <<- avg_PP_tracePlot
#                              #       rm(avg_PP_tracePlot)
                                
                
#                             #{if(x_axis_switch==TRUE)easy_remove_x_axis()}
#                             save_plot(filename=paste0("AVG_PP_waveforms_v3_hlines"),plot=avg_PP_tracePlot, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)



#### GET SPECIFIC PAIRED PULSE ISI AND ROI FOR DISPLAY OF SUBTRACTION METHOD

check_PP_df <- tmp_df %>% dplyr::filter(Ca_mM == 1, protocol == "PP60", trackROI_key == "GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009",ROI_key == "GluSnFR3-SyPhy-marker-195-dish06-plate01-region01-1Ca-PP60-repl03-ROI0009") %>% mutate(facet="1")

check_PPR_calc_peak1<- PPR_calc %>% dplyr::filter(Ca_mM == 1, protocol == "PP60", trackROI_key == "GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009",ROI_key == "GluSnFR3-SyPhy-marker-195-dish06-plate01-region01-1Ca-PP60-repl03-ROI0009") %>% mutate(facet="1")#,label=expression("Pulse")
check_PPR_calc_peak2<- PPR_calc %>% dplyr::filter(Ca_mM == 1, protocol == "PP60", trackROI_key == "GluSnFR3-SyPhy-marker-dish06-plate01-region01-ROI0009",ROI_key == "GluSnFR3-SyPhy-marker-195-dish06-plate01-region01-1Ca-PP60-repl03-ROI0009") %>% mutate(facet="3")

check_subtract_df<- check_PP_df %>% mutate(facet="3")
check_AP_df<- check_PP_df %>% mutate(facet="2")

vline_OG<- data.frame(x_int=c(0,0.06),facet=c("1","1"))
vline_AP<- data.frame(x_int=c(0),facet=c("2"))
vline_PP<-data.frame(x_int=c(0,0.06),facet=c("3","3"))
facet.labs <- c("Original Trace", "Single Stimulus Response","Subtracted Trace")
names(facet.labs) <- c(1,2,3)

facet_labels<- c("1" = "Original Trace", "2"="Test Pulse","3"="Subtracted Trace")

 subtraction_PP_plot<-ggplot(check_PP_df, aes(x=normTime, y = dFF,  colour=Ca))+
                             geom_vline(data=vline_OG,aes(xintercept=x_int),lty='dashed',colour='grey21',size=0.5)+
                             geom_vline(data=vline_AP,aes(xintercept=x_int),lty='dashed',colour='grey21',size=0.5)+
                             geom_vline(data=vline_PP,aes(xintercept=x_int),lty='dashed',colour='grey21',size=0.5)+
                             
                             geom_path(size=1,alpha=1)+
                             geom_path(data=check_AP_df,aes(x=normTime,y=singleAP_dFF), colour="grey21",size=1)+
                             geom_path(data=check_subtract_df,aes(x=normTime,y=dFF_subtracted,colour=Ca),size=1)+
                             geom_point(data=check_PPR_calc_peak1,aes(x=time_of_peak1,y=peak1_fix), shape=23, colour="blue",fill="blue",size=4,alpha=0.6)+
                             #geom_text(data=check_PPR_calc_peak1, aes(x=time_of_peak1+0.15,y=peak1_fix)) #peak1
                             geom_point(data=check_PPR_calc_peak2,aes(x=time_of_peak2,y=peak2_fix), shape=23, colour="red",fill="red",size=4,alpha=0.6)+ #peak1
                             
                             labs(#subtitle="Subtraction Routine",
                                   x="Normalized Time (s)",
                                  y=expression(Delta*"F/F"),
                                  tag="B")+#tag_vars[1])+
                             {if(color_switch)scale_colour_manual(labels=Ca.labs, values=color_override)}+
                             coord_cartesian(ylim=c(-0.25,2.5),xlim=c(-0.2,0.60))+
                             scale_x_continuous(breaks=c(0,0.5))+

                            scale_y_continuous(breaks=c(0,1,2))+
                            facet_grid(~facet, labeller = labeller(facet = function(labels)     str_wrap(facet_labels[labels], width = 10))) +
                             theme_tufte()+
                             my.theme+
                             theme(legend.position = "none")
                         #    stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.9)+
                            
                            
                         # labs( x="Normalized Time (s)",
                         #          y=expression(Delta*"F/F"),
                         #          tag="A")+#tag_vars[1])+
                         #    {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                         #    {if(color_switch)scale_colour_manual(labels=Ca.labs, values=color_override)}+
                         #    scale_x_continuous(breaks=c(0,0.5),limits=c(-0.80,0.80))+

                         #    scale_y_continuous(breaks=c(0,1,2))+
                         #    theme_tufte()+
                         #    #guides(colour="none",fill="none")+
                         #    coord_cartesian(ylim=c(0,2.5),xlim=c(-0.2,0.80))+
                         #    my.theme+
                         #    facet_grid(~protocol_expr, labeller=label_parsed)+#labeller = labeller(protocol = protocol.labs, type=label_parsed) )+ 
                         #                                    #labeller(protocol = protocol.labs, Ca = as_labeller(Ca.labs,label_parsed))+
                         #    theme(legend.position = "right",
                         #                          strip.text.y.right = element_text(angle = 0))

                            subtract_example  <<- subtraction_PP_plot
                            #        rm(avg_PP_tracePlot)
                                
                
                            save_plot(filename=paste0("draft_subtract_example"),plot=subtraction_PP_plot, plot_width=plot_height*2.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                            rm(avg_PP_tracePlot)
                            



 avg_PP_tracePlot<-ggplot(tmp_df, aes(x=normTime, y = dFF,  colour=Ca))+
                            geom_vline(xintercept = 0, lty='dashed',colour='grey21',size=0.5)+
                            geom_vline(data=vline_df,aes(xintercept = stim2_vline, group=Ca), lty='dashed',colour='grey21',size=0.5)+
                            #geom_hline(yintercept = 1, lty='dashed',colour='red',size=0.5)+
                            #geom_path(colour="grey61",size=1,alpha=0.2)+
                            stat_summary_bin(aes(group=Ca), geom="line", fun=mean, binwidth=interFrame.var,size=1.2,alpha=0.9)+
                            
                            
                         labs( x="Normalized Time (s)",
                                  y=expression(Delta*"F/F"),
                                  tag="A")+#tag_vars[1])+
                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                            {if(color_switch)scale_colour_manual(labels=Ca.labs, values=color_override)}+
                            scale_x_continuous(breaks=c(0,0.5),limits=c(-0.80,0.80))+

                            scale_y_continuous(breaks=c(0,1,2))+
                            theme_tufte()+
                            #guides(colour="none",fill="none")+
                            coord_cartesian(ylim=c(0,2.5),xlim=c(-0.2,0.80))+
                            my.theme+
                            facet_grid(~protocol_expr, labeller=label_parsed)+#labeller = labeller(protocol = protocol.labs, type=label_parsed) )+ 
                                                            #labeller(protocol = protocol.labs, Ca = as_labeller(Ca.labs,label_parsed))+
                            theme(legend.position = "right",
                                                  strip.text.y.right = element_text(angle = 0))

                            avg_peaks  <<- avg_PP_tracePlot
                            #        rm(avg_PP_tracePlot)
                                
                
                            save_plot(filename=paste0("AVG_PP_waveforms_plot"),plot=avg_PP_tracePlot, plot_width=plot_height*2.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)
                            rm(avg_PP_tracePlot)
                            



PPR_summary$Ca_expr <- factor(PPR_summary$Ca_expr, levels = c("1Ca","2Ca")) #"0pt5Ca",
levels(PPR_summary$Ca_expr) <- Ca.labs



##### FULL PPRs #######

          PPRplot<-ggplot(PPR_summary, aes(x=pairedPulse_ISI_ms, y=mean_PPR_perROI,colour=Ca))+
                                             geom_hline(yintercept=1.0, lty="dashed", linewidth=1,colour="black",alpha=0.5)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="line", fun.y = mean,size=1,alpha=1)+
                                             #geom_smooth(aes(group=Ca,colour=Ca),formula=y~x, method=y~exp(x),alpha=1,se=FALSE,size=3)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="errorbar", fun.data=mean_se,width=20,size=1, alpha=1)+
                                             stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=3,alpha=1)+
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    #title="",
                                                    tag = "C"#tag_vars[2]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(values=color_override)}+
                                            coord_cartesian(ylim=c(0.5,1.6),xlim=c(50,550))+
                                            scale_y_continuous(breaks=c(0.5,1,1.5,2.0))+
                                            scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)))+ #, guide=guide_axis(angle=45))+
                            
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position = "none")
                                    

                                    save_plot(filename=paste0("PPR_plot_v1"),plot=PPRplot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                    avg_PPR  <<- PPRplot
                                    rm(PPRplot)
                                
                
                       

##### GET PPR VIOLINS

                            PPR_summary$protocol <- factor(PPR_summary$protocol, levels = c("singleAP", "PP60", "PP75", "PP100","PP150","PP500"))
    
                                   PPR_violin<-ggplot(PPR_summary, aes(x=protocol,y=mean_PPR_perROI,colour=Ca, fill=Ca))+
                                             geom_hline(yintercept=1.0, lty="dashed", linewidth=1,colour="black",alpha=0.5)+
                                            geom_sina(size=0.5,alpha=0.7)+
                                            geom_violin(fill=NA, size=0.5,draw_quantiles = c(0.25, 0.5, 0.75),colour="black")+
                                            
                                            
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    #title="",
                                                    tag="C"#tag_vars[3]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_fill_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_fill_manual(values=color_override)}+
                                            {if(color_switch)scale_colour_manual(values=color_override,guide='none')}+
                                            coord_cartesian(ylim=c(0,4.5))+
                                            scale_y_continuous(expand=c(0,0),breaks=c(0,1,2,3,4))+
                                            scale_x_discrete(labels=c("PP60" = "60", "PP75" = "75", "PP100" = "100","PP150" = "150","PP500" = "500" ))+#, guide=guide_axis(angle=45))+
                                            #scale_x_continuous(expand=c(0,0))+
                                            
                                            facet_grid(~Ca_expr, labeller = label_parsed)+
                                            theme_tufte()+
                                            my.theme+
                                            guides(colour = guide_legend(show = FALSE),
                                                    fill = guide_legend(show = FALSE) )+
                                            theme(strip.text.y = element_blank(),
                                                    #axis.text.x = element_text(size=36),
                                                    legend.position="none")#,
                                                   #strip.text.x=element_blank())
                                    
                                     save_plot(filename=paste0("PPR_violins_v1"),plot=PPR_violin, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                    violin_PPR  <<- PPR_violin

                                    rm(PPR_violin)




# ##### GET ALL THE SPATIAL DATA ##### 



###### for loop for plotting spatial data
        ##png_dir = "Y:\\Sam/paper1_datasets/titrateCa_PP_v3/imageFiles/zstack/lowpwr_zstack/zstackOutputs_SQ"
        prefix_str = "^MAX_GluSnFR3_SyPhy_"
        suffix_str = "_0pt5Ca_zstack_lowpwr__FLATTENED.png"
        all_lut_str = "_0pt5Ca_zstack_lowpwr__FLATTENED_COLOR.png"
        dir_regex = "GluSnFR3(.*?)region\\d\\d_"
                                                                                                #saving selected protocols 
        subset_df<- PPRplot_data %>% dplyr::filter(vid_key %in% vid_keys_to_save,protocol %in% c("PP60","PP100"),Ca =="1Ca")
        #print(unique(subset_df$protocol))
        vid_key_vec <- unique(subset_df$vid_key)
        Ca_vec <- unique(subset_df$Ca)
        protocol_vec <- unique(subset_df$protocol[which(subset_df$protocol != "singleAP")])
        print(length(protocol_vec)) 
        vars_to_plot<- c("mean_PPR_perROI")
        save_bool=TRUE
        guide_bool = TRUE
        guide_limits = c(0.5,2.0)
        #guide_title = "PPR"
        levs <- c("singleAP","PP60","PP75","PP100","PP150","PP500")
        ordered_protocols <- protocol_vec[order(match(protocol_vec, levs))]
        print(ordered_protocols)


         print("Taking a crack at generating png plots")
        #subset_df<- releaseProb_calc %>% dplyr::filter(vid_key %in% vid_keys_to_save)

        ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d") 
        vid_key_vec <- unique(subset_df$vid_key)
        #vars_to_plot<- c('releaseProbability_perROI', "CV_amplitude",'mean_quanta')#c("releaseProbability_perROI")
        #save_bool=TRUE
        #guide_bool = TRUE

        #guide_limits = c(0.5,2.5)#list(c(0,1),c(0,0.5),c(0,5))
        guide_titles = c(expression(atop("PPR", "60 ms")), expression(atop("PPR","100 ms"))) #c(expression(P[iGlu]),expression(CV[Delta*"F/F"]),expression(E[Delta*"F/F"] / S[Delta*"F/F"]))# # expression(CV[tau*"decay"])#
        


        ###scale bar params! set the left-corner 
        user_x = 20
        user_y = 24.5




       #tmp_tag_var_override = ""
       tmp_tag_var_override = c("F","G") 
        print("getting ROInames") 
        ROInames <- str_extract(save_ROIs,"ROI\\d\\d\\d\\d")   
        print(ROInames)
        for (k in 1:length(vid_key_vec)){
            
                vid_df<- PPR_summary %>% dplyr::filter(vid_key == vid_key_vec[k])
                #max_gradient_val = max(vid_df$PPR)
                vid_df_for_PPR_bias<- PPR_bias_for_spatial_plot %>% dplyr::filter(vid_key == vid_key_vec[k])

                print(vid_df_for_PPR_bias)
             

                count=0
                for(i in 1:length(Ca_vec)) {
                                    
                                    tmp_plot = NULL
                                    Ca_df<- vid_df %>% dplyr::filter(Ca == Ca_vec[i])
                                    Ca_df_for_PPR_bias<- vid_df_for_PPR_bias %>% dplyr::filter(Ca == Ca_vec[i], protocol != "singleAP")
                                    special_guide_limits = c(60,500)
                                    special_guide_titles = c("PPR bias")#c(expression(P[iGlu]),expression(CV[Delta*"F/F"]),expression(E[Delta*"F/F"] / S[Delta*"F/F"]))# # expression(CV[tau*"decay"])#
                                    special_var_to_plot = c("PPR_bias_numeric")
                                    special_tag = c("E")#c("E")
                                    special_scalebar = TRUE
                                    special_guide_bool = TRUE

                                    
                                    print("I hope this works! PNG plot of PPR bias")

                                     tmp_plot = png_plotter(df=Ca_df_for_PPR_bias,
                                                                    png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                                                    vars_to_plot=special_var_to_plot, 
                                                                    save_bool=save_bool,guide_bool=special_guide_bool,tag_var =special_tag,
                                                                    ROI_IDs = ROInames,guide_title=special_guide_titles,guide_lim=special_guide_limits,user_x = user_x, user_y = user_y,scalebar_switch=special_scalebar,add_arrows=TRUE)

                                    check_data<- Ca_df %>% select(Ca, protocol, mean_PPR_perROI)
                                    #print(check_data)
                                            
                                                for (j in 1:length(ordered_protocols)) {
                                                    protocol_df <- Ca_df %>% dplyr::filter(protocol == ordered_protocols[j])
                                                    #protocol_df_for_PPR_bias<- vid_df_for_PPR_bias %>% dplyr::filter(Ca == Ca_vec[i])

                                                            count = count+1
                                                            print(paste0(count, " was value of iteration counter at: Ca = ", unique(Ca_df$Ca), " and protocol = ", unique(protocol_df$protocol) ) )
                                                            print("This is the current dataset:")
                                                            print(check_data)
                                                                    if(save_bool==TRUE) { 
                                                                        guide_bool = TRUE
                                                                        var_to_plot = vars_to_plot
                                                                        tmp_guide_title = guide_titles[j]
                                                                        tmp_limits = guide_limits
                                                                        tmp_tag = tmp_tag_var_override[j]
                                                                        if(j == 1){tmp_scalebar = FALSE#TRUE 
                                                                                    guide_bool=TRUE
                                                                                    } else {tmp_scalebar = FALSE 
                                                                                            guide_bool=TRUE
                                                                                               }
                                                            

                                        
                                           tmp_plot = png_plotter(df=protocol_df,
                                                                    png_dir=png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
                                                                    vars_to_plot=var_to_plot, 
                                                                    save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
                                                                    ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar,add_arrows=TRUE)
                                        }

                                    }

                                                                    

                }

        }
                                    


                        
print("finished plotting maps")

# ##### step 1. establish subset of the dataframe that can be used to plot traces and physiology

                        
tmp_trace_data <- tmp_df %>% dplyr::filter(protocol_expr != "Test~Pulse") #trackROI_key %in% ROIs_to_save,
                                           #protocol %in% c("PP60","PP75","PP100","PP150","PP500"))
trace_tags = c("H","J","P") #c("B","D")#
PPR_tags = c("I","K","Q") #c("C","E")#
pulse1_tags = c("J","M","R")
pulse2_tags = c("K","N","S")
tracked_ROIs <- ROIs_to_save#unique(tmp_trace_data$trackROI_key)
print(tracked_ROIs)
print(ROIs_to_save)

            for (i in 1:length(ROIs_to_save)) {
                    print(paste0("loop iter on count: ", i))
                    
                    ROI_title<-gsub("GluSnFR3-SyPhy-marker-dish\\d\\d-plate\\d\\d-region\\d\\d-","",tracked_ROIs[i])
                    tmp_summary_data <- PPRplot_data %>%
                                        dplyr::filter(trackROI_key == ROIs_to_save[i])
                    tmp_trace_subset<- tmp_trace_data %>% 
                                        dplyr::filter(trackROI_key == ROIs_to_save[i]) %>%
                                        mutate(stim2_vline = as.numeric(gsub("PP", "", protocol) )/1000 )
                    keep_cols <- c('Ca','protocol','Ca_expr','protocol_expr','stim2_vline')
                    x_ints <- tmp_trace_subset[,(names(tmp_trace_subset) %in% keep_cols)]
                    x_ints <- unique(x_ints)
                    


                    #print(unique(x_ints$stim2_vline)) 

                    #print("these should be the x-intercepts")
                    #print(x_ints)
                    #tmp_trace_subset$Ca<- as.factor(as.character( (tmp_trace_subset$Ca) ) )
                    #levels(tmp_trace_subset$Ca) <- new_levels
                        
                        

                          tmp_PPRplot<-ggplot(tmp_summary_data, aes(x=pairedPulse_ISI_ms, y=PPR,colour=Ca))+
                                             geom_hline(yintercept = 1, size=1,lty="dashed",colour="black",alpha=0.5)+
                                             stat_summary(geom="line", fun.y = mean, size=1, alpha=0.85 )+
                                             geom_point(alpha=0.5, size=2)+
                                             
                                             
                                             stat_summary(geom="errorbar", fun.data=mean_se, width=10, size=1, alpha=1)+
                                             stat_summary(geom="point", fun.y=mean, size=2.5,shape=21,fill="white",stroke=1.5,alpha=1)+
                                            #stat_summary(aes(group=Ca,colour=Ca),geom="line", fun.y=mean,alpha=1,size=1.5,na.rm=TRUE)+
                                            #stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=3,alpha=1,na.rm=TRUE)+
                                             
                                            labs( x="ISI (ms)",
                                                    y=expression(Pulse[2]/Pulse[1]),
                                                    #title="",
                                                    tag = PPR_tags[i]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual(values=color_override)}+
                                            coord_cartesian(ylim=c(0.5,2.5),xlim=c(50,500))+ 
                                            scale_y_continuous(breaks=c(0.5,1,2))+
                                            scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)))+
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position="none")

                                            save_plot(filename=paste0("PPR_plot_", tracked_ROIs[i],"_"),plot=tmp_PPRplot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                            nam <- paste("tmp_PPRplot", i, sep = "")
                                            assign(nam, tmp_PPRplot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                        tmp_pulse_plot<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=peak1_fix,colour=Ca))+
                                             geom_hline(yintercept = 0.4, size=0.7,lty="dashed",colour="black")+
                                             stat_summary(geom="line", fun.y = mean, size=1, alpha=0.85,colour='black' )+
                                             geom_point(alpha=0.5, size=2)+
                                             
                                             
                                             stat_summary(geom="errorbar", fun.data=mean_se, width=0.1, size=1, alpha=1,colour='black')+
                                             stat_summary(geom="point", fun.y=mean, size=2.5,alpha=1,colour='black')+
                                            #stat_summary(aes(group=Ca,colour=Ca),geom="line", fun.y=mean,alpha=1,size=1.5,na.rm=TRUE)+
                                            #stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=3,alpha=1,na.rm=TRUE)+
                                             
                                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                                    y=expression(Pulse[1]~Delta*"F/F"),
                                                    #title="",
                                                    tag = pulse1_tags[i]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual( values=color_override)}+
                                            coord_cartesian(ylim=c(0,4),xlim=c(0.25,2.25))+ 
                                            scale_y_continuous(breaks=c(0,2,4))+
                                            scale_x_continuous(breaks=c(0.5,1,2),labels=c("0.5","1","2"))+#scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)))+
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position="none")

                                            save_plot(filename=paste0("PPR_plot_pulse1_amp_", tracked_ROIs[i],"_"),plot=tmp_pulse_plot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                            nam <- paste("tmp_pulse1_plot", i, sep = "")
                                            #assign(nam, tmp_pulse_plot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                                            #nam <- paste("tmp_PPRplot", i, sep = "")
                                            #assign(nam, tmp_PPRplot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                         tmp_pulse_plot<-ggplot(tmp_summary_data, aes(x=Ca_mM, y=peak2_fix,colour=Ca))+
                                             geom_hline(yintercept = 0.4, size=0.7,lty="dashed",colour="black")+
                                             stat_summary(geom="line", fun.y = mean, size=1, alpha=0.85,colour='black' )+
                                             geom_point(alpha=0.5, size=2)+
                                             
                                             
                                             stat_summary(geom="errorbar", fun.data=mean_se, width=0.1, size=1, alpha=1,colour='black')+
                                             stat_summary(geom="point", fun.y=mean, size=2.5,alpha=1,colour='black')+
                                            #stat_summary(aes(group=Ca,colour=Ca),geom="line", fun.y=mean,alpha=1,size=1.5,na.rm=TRUE)+
                                            #stat_summary(aes(group=Ca, colour=Ca), geom="point", fun.y=mean,size=3,alpha=1,na.rm=TRUE)+
                                             
                                            labs( x=expression("["*Ca^'2+'*"] (mM)"),
                                                    y=expression(Pulse[2]~Delta*"F/F"),
                                                    #title="",
                                                    tag = pulse2_tags[i]
                                                    )+
                                            
                                            {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                            {if(color_switch)scale_colour_manual( values=color_override)}+
                                            coord_cartesian(ylim=c(0,4),xlim=c(0.25,2.25))+ 
                                            scale_y_continuous(breaks=c(0,2,4))+
                                            scale_x_continuous(breaks=c(0.5,1,2),labels=c("0.5","1","2"))+#scale_x_continuous(breaks=c(60,75,100,150,500), labels = remove_labels(scales::label_number_auto(), c(60,100)))+
                                            theme_tufte()+
                                            my.theme+
                                            theme(legend.position="none")

                                            save_plot(filename=paste0("PPR_plot_pulse2_amp_", tracked_ROIs[i],"_"),plot=tmp_pulse_plot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                            nam <- paste("tmp_pulse2_plot", i, sep = "")
                                            #assign(nam, tmp_pulse_plot,envir = .GlobalEnv)#spat_facet <- tmp_plot





                        # tmp_trace_avgs<-tmp_trace_subset %>% ungroup() %>% select(Ca,protocol,reindex, avg_normTime,avg_dFF) %>% group_by(Ca,protocol) %>% arrange(reindex)
                        # tmp_trace_avgs<-unique(tmp_trace_avgs)
                        # print(tmp_trace_avgs)            
                        print(names(tmp_trace_subset))        

                        tmp_tracePlot<- ggplot(tmp_trace_subset, aes(x=normTime, y = dFF,  group=exposeNum, colour=Ca,alpha=Ca))+
                                                    geom_vline(xintercept = 0, lty="dashed",colour='grey21',size=0.75)+
                                                    geom_vline(data=x_ints, aes(xintercept = stim2_vline,group=Ca), lty="dashed",colour='grey21',size=0.75)+
                                                    
                                                    geom_path(data=tmp_trace_subset,aes(x=normTime,y=dFF,group=exposeNum,colour=Ca,alpha=Ca),size=1)+
                                                    stat_summary_bin(aes(group=Ca), colour="black",geom="line", fun=mean, binwidth=interFrame.var,size=0.5,alpha=0.95)+
                                                    
                                                    #geom_path(data=tmp_trace_avgs, aes(x=avg_normTime,y=avg_dFF,group=interaction(Ca,protocol)), colour="black",size=0.9,alpha=0.95)+
                                                    
                                                    labs( x="Normalized time (s)",
                                                          y=expression(Delta*"F/F"),
                                                          title=paste0("Bouton",gsub("ROI000", "  ", ROI_title)),
                                                          #colour=expression("[Ca"^{"2+"}*"], mM"),
                                                                        tag = trace_tags[i])+
                                                    {if(color_switch == FALSE)scale_colour_brewer(palette = "Dark2")}+
                                                    {if(color_switch)scale_colour_manual(values=color_override)}+
                                                    {if(color_switch)scale_alpha_manual(values=c(0.8,0.7,0.7))}+
                                                    scale_y_continuous(breaks=c(0,2,4,6))+
                                                    scale_x_continuous(breaks=c(0,0.5))+
                                                    theme_tufte()+
                                                    guides(colour="none",fill="none")+
                                                    coord_cartesian(ylim=c(-0.25,4.5),xlim=c(-0.2,0.80))+
                                                    my.theme+
                                                    facet_grid(Ca_expr~protocol_expr, labeller = label_parsed)+# labeller(protocol = protocol.labs, Ca = Ca.labs, type=label_parsed))+
                                                    theme(legend.position = "none",
                                                            strip.text.y.right = element_blank()
                                                            #strip.text.x = element_text(colour="black",size=22*scalefactor,family="sans")
                                                            )

                                  #  save_plot(filename=paste0("PairedPulse_traces_", tracked_ROIs[i],"_"),plot=tmp_tracePlot, plot_width=plot_height,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)

                                    nam <- paste("tmp_tracePlot", i, sep = "")
                                    assign(nam, tmp_tracePlot,envir = .GlobalEnv)#spat_facet <- tmp_plot

                                            

                        # gs = list(tmp_tracePlot,  #1
                        #          tmp_PPRplot)  #3
                        # margin = theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))              
                    



                        # hlay <- rbind( c(1,1,1,1,1,1,1,NA,NA,NA,NA,NA),
                        #                c(1,1,1,1,1,1,1,2,2,2,2,2),
                        #                c(1,1,1,1,1,1,1,2,2,2,2,2),
                        #                c(1,1,1,1,1,1,1,2,2,2,2,2),
                        #                c(1,1,1,1,1,1,1,2,2,2,2,2),
                        #                c(1,1,1,1,1,1,1,2,2,2,2,2),
                        #                c(1,1,1,1,1,1,1,NA,NA,NA,NA,NA)
                        #            )
                
                    

                    
                        # ROI_report = grid.arrange(grobs=lapply(gs, "+", margin), layout_matrix=hlay)
                        # #ggsave(filename=paste0(tracked_ROIs[i],"_releaseProb.emf"),plot=ROI_report,dpi=900, units="in",width=24,height=18, device = {function(filename, ...) devEMF::emf(file = filename, ...)})
                        # save_plot(filename=paste0("ROI_report_", tracked_ROIs[i],"_"),plot=ROI_report, plot_width=plot_height*1.5,plot_height=plot_height, scale_f = scalefactor,dpi_val = 600)



    }



  # tmp_trace_avgs
 
 output_data<- list(PPR_calc_all_obs = PPR_removedup, PPR_perROI = PPR_summary, PPR_avgs = plot_avgs,pie_data = pie_chart_final,PPR_bias = PPR_bias_final, PPR_bias_spatial_plot = PPR_bias_for_spatial_plot,tmp_trace_subset = tmp_trace_subset, ROIs_to_save = ROIs_to_save, tmp_summary = tmp_summary_data, tmp_df = tmp_df )
 output_data
#PPR_bias_for_spatial_plot
#check_PP_df

}						
						