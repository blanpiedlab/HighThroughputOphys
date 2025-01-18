##transmission_stats_v2.R
## Written by Samuel T. Barlow
## 04.17.23

# transmission_stats_v2.R is a module designed to extract meaningful metrics from dual-color synaptic transmission imaging data. 

# 1. Synaptic transmission rate as a function of chemical_condition (APV vs. ctrl)
# 2. JF646 and GluSnFR3 spike properties (by groups chemical_condition, as well as subsetted into failures and successes)
# 3. generate 2D plot of GluSnFR3 spikes vs. JF646amplitude
# 4. generate averaged traces which exploit a windowed view of 


subDir <- "transmission_stats_t2"   


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


transmission_stats<- function(traces_df, interSpike_thresh = NULL, peak_groupers,groupers, positive_ROI_groupers,plotBy, levels, secondAxis = NULL, color_override=color_override,comparisons, save_ROIs = ROIs_to_save){

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

##### first, generate a dataframe that achieves the following: 
# A. Look at tracked ROIs
# B. Begin in the control condition
# C. Define ctrl as "1" and vehicle or APV as "2"
# D. Filter out all ROIs that have <1 event in ctrl, be that Glu or JF646. 
                
                traces_df_fix_brokenROI<- traces_df %>% dplyr::filter(chemical_condition == "APV", vid_key == "dish7-plate07-region04", absoluteTime >= 6) ###this dataset had a focus issue that is messing with quantification
                traces_df_revise_brokenROI <- traces_df %>% dplyr::filter(chemical_condition == "ctrl",vid_key == "dish7-plate07-region04")
                traces_df_omit_brokenROI<- traces_df %>% dplyr::filter(vid_key != "dish7-plate07-region04" )
                traces_df_fixed <- bind_rows(traces_df_omit_brokenROI, traces_df_fix_brokenROI,traces_df_revise_brokenROI)


                variance_filter<- traces_df_fixed %>% dplyr::filter(alpha_str == "notPeak",chemical_condition == "ctrl",sensor == "JF646") %>% 
                                                      group_by(trackROI_key) %>%
                                                      summarise(sd_trace = sd(dFF,na.rm=TRUE)) %>%
                                                      ungroup() %>%
                                                      dplyr::filter(sd_trace >= quantile(sd_trace, 0.85, na.rm=TRUE))
                traces_df_fixed<- traces_df_fixed %>% dplyr::filter(!trackROI_key %in% unique(variance_filter$trackROI_key) )


                get_peak_maxima <- traces_df_fixed %>% 
                                             group_by_at(peak_groupers) %>% 
                                             summarise(variance=var(dFF,na.rm=TRUE),
                                                        noise = sd(dFF, na.rm=TRUE),
                                                        mean_noise = mean(dFF[which(peakID == "NotPeak")], na.rm=TRUE),
                                                        signal = max(dFF) ) %>%
                                            ungroup() %>%
                                            group_by_at(positive_ROI_groupers) %>%
                                            summarise(variance=variance[which(peakID=="NotPeak")],
                                                       noise = abs(noise[which(peakID == "NotPeak")]),
                                                       mean_noise=mean_noise[which(peakID=="NotPeak")]#,
                                                       #min_signal =ifelse(is.finite(signal[which(peakID !="NotPeak")]),min(signal[which(peakID !="NotPeak")],na.rm=TRUE), NA),
                                                       #max_signal =ifelse(is.finite(signal[which(peakID !="NotPeak")]),max(signal[which(peakID !="NotPeak")],na.rm=TRUE), NA),
                                                       #min_SNR = min_signal/noise,
                                                       #SNRpass = ifelse(is.na(noise), FALSE, min_SNR >= 2) #ridiculous noise to filter on but we ball
                                                       ) %>%
                                                       dplyr::filter(variance < 0.02)
                just_maxima <- traces_df_fixed %>% group_by_at(peak_groupers) %>%
                                             dplyr::filter(peakID != "NotPeak") %>% 
                                             summarise(amplitude = max(dFF) )
                
                just_noise<-traces_df_fixed %>% group_by_at(peak_groupers) %>% dplyr::filter(peakID=="NotPeak") %>% rename(noise = dFF) %>% mutate(groupID = paste(sensor,chemical_condition,sep='-'))
                                           


                 get_peak_sums <- suppressMessages(traces_df_fixed %>% 
                                                         group_by_at(positive_ROI_groupers) %>%
                                                         #dplyr::filter(peakID != "NotPeak") %>% 
                                                         summarise(peaklist= paste(unique(peakID), collapse="-"), #doesn't seem to be throwing NAs, so trackROIs are getting filtered when they have no peaks identified. working as intended.
                                                                    peak_sum = length(unique(peakID))-1,
                                                                    time_dur = max(absoluteTime, na.rm=TRUE) - min(absoluteTime,na.rm=TRUE),
                                                                    peak_freq = peak_sum/time_dur) %>%
                                                         mutate(exchange_status = case_when(chemical_condition == "ctrl" ~ "Before",
                                                                                            chemical_condition == "vehicle" | chemical_condition == "APV" ~ "After"),
                                                                    activity_ID = case_when( #( peak_sum >=15 & sensor == "GluSnFR3")  ~ "Glu_overactive",
                                                                                            #( peak_sum >=15 & sensor == "JF646")  ~ "JF_overactive",
                                                                                                ( peak_sum >=1 & sensor == "GluSnFR3")  ~ "Glu_active",
                                                                                                ( peak_sum == 0 & sensor == "GluSnFR3") ~ "Glu_inactive",
                                                                                                ( peak_sum >=1 & sensor == "JF646")  ~ "JF_active",
                                                                                                (peak_sum == 0 & sensor == "JF646") ~ "JF_inactive") ) #,
                                                                                                #(chemical_condition == "vehicle" | chemical_condition == "APV") ~ "exchange") )  
                                                         )

                cleanup_ROIs <- suppressMessages(left_join(get_peak_sums,get_peak_maxima) %>%
                                                            ungroup() %>%  
                                                            group_by(trackROI_key) %>%
                                                            dplyr::filter(exchange_status == "After",activity_ID == "Glu_active") #%>%SNRpass == TRUE,
                                                            #mutate(exchange_label=paste(sort(unique(exchange_status)), collapse="&")) %>% dplyr::filter(sensorIDs == "Before&After")

                                                             )

                clean_ROI_list <- unique(cleanup_ROIs$trackROI_key)
                

#### having generated a summary dataframe that keeps all of the sensor_trackROI_keys that feature at least 1 event, we can now do estimates of frequency on the basis of chemical_condition (control or vehicle/APV)
            
                 clean_freq <- get_peak_sums %>% 
                                            ungroup() %>%
                                            group_by(sensor_trackROI_key) %>%
                                            dplyr::filter(trackROI_key %in% clean_ROI_list) %>% 
                                            mutate(condition_ID = case_when(chemical_condition == "ctrl" ~ 1,
                                                                            chemical_condition == "APV" ~ 2,
                                                                            chemical_condition == "vehicle" ~2,
                                                                            chemical_condition == "100serine" ~ 2,
                                                                            chemical_condition == "100kyna" ~ 2),
                                                    baseline_freq = peak_freq[which(condition_ID == 1)],
                                                    norm_freq = peak_freq/baseline_freq
                                                    )

                clean_vehicle_ROIs<- clean_freq %>% dplyr::filter(chemical_condition == "vehicle") 
                clean_vehicle_ROI_list<- unique(clean_vehicle_ROIs$trackROI_key)
                clean_freq_vehicle <- clean_freq %>% dplyr::filter(trackROI_key %in% clean_vehicle_ROI_list) %>% mutate(groupID = paste(sensor,"Vehicle",sep='-'),facetID = "Control" )
                
                clean_APV_ROIs<- clean_freq %>% dplyr::filter(chemical_condition == "APV") 
                clean_APV_ROI_list<- unique(clean_APV_ROIs$trackROI_key)
                clean_freq_APV <- clean_freq %>% dplyr::filter(trackROI_key %in% clean_APV_ROI_list) %>% mutate(groupID = paste(sensor,"AP5",sep='-'),facetID = "AP5" )

                clean_freq_revised <- bind_rows(clean_freq_vehicle, clean_freq_APV) %>% ungroup() %>% group_by(trackROI_key) %>% mutate(sensorIDs=paste(sort(unique(sensor)), collapse="&")) %>% dplyr::filter(sensorIDs == "GluSnFR3&JF646")
                #


                APV_labels = c("Control", "AP5")
                vehicle_labels = c('Control', "Veh")

# ### GENERATE A PLOT OF PEAK FREQUENCY AS A FUNCTION OF EXCHANGE

                    
                    levels_groupID = c('GluSnFR3-Vehicle','JF646-Vehicle','GluSnFR3-AP5','JF646-AP5')
                    clean_freq_revised$groupID=factor(clean_freq_revised$groupID, levels= levels_groupID)
                    new_levels_groupID = c('iGluSnFR3-Control','JF646-Control','iGluSnFR3-AP5','JF646-AP5')
                    levels(clean_freq_revised$groupID) = new_levels_groupID
                    clean_freq_revised$facetID=factor(clean_freq_revised$facetID, levels=c("Control","AP5"),labels=c("Control", expression("100"~mu*"M AP5")))
                    #levels(clean_freq_revised$facetID
                    labels_revised = c('Before', 'After')
                    peak_sum_plot_v2 <- ggplot(clean_freq_revised, aes(x =  condition_ID, y = peak_freq, alpha=groupID) ) +
                                           #geom_line(aes(group =trackROI_key), alpha=0.05, colour="grey21", linewidth=1)+
                                           #geom_point(aes(group = trackROI_key), alpha=0.05, colour="grey21", size=2)+
                       

                                           stat_summary(aes(group=groupID, colour=groupID,linetype=groupID), geom="line", fun=mean,linewidth=2,alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID), geom="errorbar", fun.data=mean_se,width=0.1,size=2, alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID,shape=groupID), geom="point", stroke=2,fill='white',fun=mean,size=4,alpha=1,na.rm=TRUE)+
                                           labs( x="",
                                                  y=expression("Peak frequency ("*s^{-1}*")"),
                                                  #title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="H")+
                                           scale_colour_manual(name = "", values = c("iGluSnFR3-Control" = "#339900", 'iGluSnFR3-AP5' = '#339900', "JF646-Control" ="#FF33FF","JF646-AP5" ="#FF33FF"))+
                                           scale_shape_manual(name = "", values = c("iGluSnFR3-Control" = 16, 'iGluSnFR3-AP5' = 21, "JF646-Control" =16,"JF646-AP5" = 21))+
                                           scale_linetype_manual(name = "", values = c("iGluSnFR3-Control" = "solid", 'iGluSnFR3-AP5' = 'longdash', "JF646-Control" ="solid","JF646-AP5" ="longdash"))+
                                           scale_alpha_manual(name="",values = c("iGluSnFR3-Control" = 1.0, "iGluSnFR3-AP5" = 0.8, "JF646-Control" = 1.0, "JF646-AP5" =0.8))+
                                           
                                            coord_cartesian(ylim=c(0.05,0.35))+#xlim=c(0.9,2.1),
                                           #scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
                                           scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~facetID,labeller=label_parsed)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_text(colour="black", size=24, family="sans"),
                                                 legend.position = c(0.20,0.15),#"right",
                                                 legend.justification="center"
                                            )



#### save to global env
                                peak_rate_H  <<- peak_sum_plot_v2

         
                            ggsave(filename=paste0("peak_sum_plot_v2.png"),plot=peak_sum_plot_v2, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
        

###### GENERATE A VIOLIN PLOT OF PEAK AMPLITUDE AS A FUNCTION OF EXCHANGE

                        
                    


                          clean_maxima_APV <- just_maxima %>% dplyr::filter(trackROI_key %in% clean_APV_ROI_list) %>%
                                                                mutate(condition_ID = case_when(chemical_condition == "ctrl" ~ 1,
                                                                                                chemical_condition == "APV" ~ 2,
                                                                                                chemical_condition == "vehicle" ~2,
                                                                                                chemical_condition == "100serine" ~ 2,
                                                                                                chemical_condition == "100kyna" ~ 2),
                                                                        groupID = paste(sensor,chemical_condition,sep='-')
                                                                )

                   


                          clean_maxima_vehicle <- just_maxima %>% dplyr::filter(trackROI_key %in% clean_vehicle_ROI_list) %>%
                                                                mutate(condition_ID = case_when(chemical_condition == "ctrl" ~ 1,
                                                                                                chemical_condition == "APV" ~ 2,
                                                                                                chemical_condition == "vehicle" ~2,
                                                                                                chemical_condition == "100serine" ~ 2,
                                                                                                chemical_condition == "100kyna" ~ 2),
                                                                        groupID = paste(sensor,chemical_condition,sep='-')
                                                                )

                          clean_maxima_revised <- bind_rows(clean_maxima_vehicle,clean_maxima_APV) %>%  
                                                                ungroup() %>%
                                                                group_by(sensor,chemical_condition,groupID) %>%
                                                                dplyr::filter(amplitude>0.03) %>% 
                                                                mutate(med_amplitude = median(amplitude,na.rm=TRUE) ) 
                            
                          clean_noise<- just_noise %>% dplyr::filter(trackROI_key %in% clean_ROI_list) 
                          

                          levels= c("GluSnFR3-ctrl", "GluSnFR3-vehicle", 'GluSnFR3-APV',"JF646-ctrl", "JF646-vehicle","JF646-APV" )
                          clean_maxima_revised$groupID = factor(clean_maxima_revised$groupID, levels= levels)
                          chemical_levels = c('ctrl','vehicle','APV')
                          clean_maxima_revised$chemical_condition = factor(clean_maxima_revised$chemical_condition, levels= chemical_levels, labels = c("Pre-Treatment", "Vehicle", "AP5"))
                          clean_noise$chemical_condition = factor(clean_noise$chemical_condition, levels= chemical_levels, labels = c("Pre-Treatment", "Vehicle", "AP5"))
                          
                          lower_bound_x = 0
                          upper_bound_x = 2.5
                          binwidth = (upper_bound_x - lower_bound_x) / 63  

                          amplitude_histogram<-ggplot(clean_maxima_revised, aes(x=amplitude, group=groupID,colour=groupID,fill=groupID,alpha=groupID))+
                                                #geom_histogram(data=clean_noise, aes(x=noise, y=after_stat(ncount),group=groupID),fill='grey21',binwidth = binwidth, alpha=0.5,boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                #geom_freqpoly(data=clean_noise, aes(x=noise,y=after_stat(ncount)),colour='grey21',binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5)+
                                                
                                                geom_histogram(aes(y=after_stat(ncount)),binwidth = binwidth, boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                geom_freqpoly(aes(y=after_stat(ncount)),binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5)+
                                                
                                                geom_vline(aes(xintercept = med_amplitude), colour= "black", lty = "dashed", linewidth = 2,alpha=1)+
                                      labs( x=expression("Signal amplitude ("*Delta*"F/F)"),
                                                  y=expression(N/N[max]),
                                                  #title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="D")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("GluSnFR3-vehicle" = "#336600", "GluSnFR3-APV" = "#669900", 'GluSnFR3-ctrl' = '#339900',"JF646-vehicle" = "#990066", "JF646-APV" ="#993366","JF646-ctrl" ="#FF33FF"))+
                                           scale_fill_manual(name = "", values = c("GluSnFR3-vehicle" = "#336600", "GluSnFR3-APV" = "#669900", 'GluSnFR3-ctrl' = '#339900',"JF646-vehicle" = "#990066", "JF646-APV" ="#993366","JF646-ctrl" ="#FF33FF"))+
                                           scale_alpha_manual(name="",values = c("GluSnFR3-vehicle" = 0.5, "GluSnFR3-APV" = 0.4, 'GluSnFR3-ctrl' = 0.8,"JF646-vehicle" = 0.5, "JF646-APV" =0.4,"JF646-ctrl" =0.8))+
                                           coord_cartesian(ylim=c(0,1.1),xlim=c(0,2.0))+
                                           scale_y_continuous(expand=c(0,0))+
                                           facet_grid(chemical_condition~sensor,scales="free_x")+
                                           scale_x_continuous(limits=c(0,2.5),breaks=c(0,0.5,1.0,1.5,2.0,2.5))+#scale_y_continuous(limits=c())
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank(),
                                                 strip.text.y=element_text(angle = 0, size=28, family="sans"),
                                                 panel.spacing = unit(1, "cm"),
                                                 legend.position="none"
                                            )
                                

#### save to global env
                                #histo_G  <<- amplitude_histogram


                               save_plot(filename=paste0("histo_plot_amplitude_v2"),plot=amplitude_histogram, plot_width=plot_dim*1.5,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

                                amplitude_histogram<-ggplot(clean_maxima_revised, aes(x=amplitude, group=groupID,colour=groupID,fill=groupID,alpha=groupID))+
                                                geom_histogram(data=clean_noise, aes(x=noise, y=after_stat(ncount),group=groupID),fill='grey21',binwidth = binwidth, alpha=0.5,boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                geom_freqpoly(data=clean_noise, aes(x=noise,y=after_stat(ncount)),colour='grey21',binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5)+
                                                
                                                geom_histogram(aes(y=after_stat(ncount)),binwidth = binwidth, boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                geom_freqpoly(aes(y=after_stat(ncount)),binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5)+
                                                
                                                geom_vline(aes(xintercept = med_amplitude), colour= "black", lty = "dashed", linewidth = 2,alpha=1)+
                                      labs( x=expression("Signal amplitude ("*Delta*"F/F)"),
                                                  y=expression(N/N[max]),
                                                  #title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="B")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("GluSnFR3-vehicle" = "#336600", "GluSnFR3-APV" = "#669900", 'GluSnFR3-ctrl' = '#339900',"JF646-vehicle" = "#990066", "JF646-APV" ="#993366","JF646-ctrl" ="#FF33FF"))+
                                           scale_fill_manual(name = "", values = c("GluSnFR3-vehicle" = "#336600", "GluSnFR3-APV" = "#669900", 'GluSnFR3-ctrl' = '#339900',"JF646-vehicle" = "#990066", "JF646-APV" ="#993366","JF646-ctrl" ="#FF33FF"))+
                                           scale_alpha_manual(name="",values = c("GluSnFR3-vehicle" = 0.5, "GluSnFR3-APV" = 0.4, 'GluSnFR3-ctrl' = 0.8,"JF646-vehicle" = 0.5, "JF646-APV" =0.4,"JF646-ctrl" =0.8))+
                                            coord_cartesian(,ylim=c(0,1.1))+
                                           scale_y_continuous(expand=c(0,0), breaks=c(0,0.5,1.0))+
                                           facet_grid(chemical_condition~sensor,scales="free_x")+
                                           scale_x_continuous(breaks=c(0,0.5,1.0,1.5,2.0,2.5))+#scale_y_continuous(limits=c())
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank(),
                                                 strip.text.y=element_text(angle = 0, size=28, family="sans"),
                                                 panel.spacing = unit(1, "cm"),
                                                 legend.position="none"
                                            )
                 
                               #ggsave(filename=paste0("histo_plot_amplitude_with_noise_v2.png"),plot=amplitude_histogram, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)


### get raw counts of peaks
                    count_peaks <- clean_freq_revised %>% group_by(sensor,chemical_condition) %>% summarise(sum_peaks = sum(peak_sum),
                                                                                            sum_ROIs = length(unique(trackROI_key)),
                                                                                            sum_videos = length(unique(vid_key) ),
                                                                                            norm_peak_sum = sum_peaks/sum_ROIs) %>%
                                                                                            mutate(groupID = paste(sensor,chemical_condition,sep='-'))
                    #count_peaks<-unique(count_peaks_v0)

                    count_peaks$groupID = factor(count_peaks$groupID, levels= levels)
                    count_peaks$chemical_condition = factor(count_peaks$chemical_condition, levels= chemical_levels, labels = c("Pre-Treatment", "Vehicle", "AP5"))

                     count_dots<-ggplot(count_peaks, aes(x=chemical_condition, y=sum_peaks,group=groupID,colour=groupID,shape=groupID,alpha=groupID))+
                                           geom_point(size=8,stroke=3,fill='white')+
                                           labs( x="",
                                                  y=expression("Total Peaks Sampled"),
                                                  #title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="B")+
                                           scale_colour_manual(name = "", values = c("GluSnFR3-vehicle" = "#336600", "GluSnFR3-APV" = "#669900", 'GluSnFR3-ctrl' = '#339900',"JF646-vehicle" = "#990066", "JF646-APV" ="#993366","JF646-ctrl" ="#FF33FF"))+
                                           scale_shape_manual(name = "", values = c("GluSnFR3-vehicle" = 16, "GluSnFR3-APV" = 21, 'GluSnFR3-ctrl' = 16,"JF646-vehicle" = 16, "JF646-APV" =21,"JF646-ctrl" =16))+
                                           scale_alpha_manual(name="",values = c("GluSnFR3-vehicle" = 0.5, "GluSnFR3-APV" = 0.4, 'GluSnFR3-ctrl' = 0.8,"JF646-vehicle" = 0.5, "JF646-APV" =0.4,"JF646-ctrl" =0.8))+
                                           
                                            #coord_cartesian(xlim=c(0.9,2.1),ylim=c(1,3.25))+
                                           #scale_y_continuous(breaks=c(0,1,2,3,4,5))+
                                           #scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~sensor)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank()
                                            )
                 
                            ggsave(filename=paste0("counts_sum_peaks_v1.png"),plot=count_dots, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
        
                           


###### GET SYNAPTIC TRANSMISSION INFO ####### 
                    
                    #mutate call here is a handy way of determining if 2 different strings are present in a single group. Identify all the ROIs (synapses) where we can measure synaptic transmission because both sensor traces exist. 
                   # cleanup_ROIs_both<- cleanup_ROIs %>% select(trackROI_key,sensor) %>% group_by(trackROI_key) %>% mutate(sensorIDs=paste(sort(unique(sensor)), collapse="&")) %>% dplyr::filter(sensorIDs == "GluSnFR3&JF646")
                    #clean_ROI_list_both<- unique(cleanup_ROIs_both$trackROI_key)

                    #alternative approach to getting ROIs that feature a lot of activity
                    #active_ROIs_Glu<-clean_freq_revised %>% dplyr::filter(peak_sum >= 1,chemical_condition == "ctrl",sensor == "GluSnFR3") 
                    #active_ROIs_Glu_list <- unique(active_ROIs_Glu$trackROI_key) 
                    #activeROIs_both <- clean_freq_revised %>% dplyr::filter(trackROI_key %in% active_ROIs_Glu_list) %>% ungroup() %>% select(trackROI_key,sensor) %>% group_by(trackROI_key) %>% mutate(sensorIDs=paste(sort(unique(sensor)), collapse="&")) %>% dplyr::filter(sensorIDs == "GluSnFR3&JF646")
                    #clean_activeROI_list_both<- unique(activeROIs_both$trackROI_key)
                    #subset_cleanROIs<- clean_freq %>% dplyr::filter(trackROI_key %in% clean_ROI_list) %>% select(trackROI_key)
                    #subset_cleanROIs_list<-unique(subset_cleanROIs$trackROI_key)
                    total_ROIs<- length(clean_ROI_list)
                    ROIs_to_sample = ceiling( total_ROIs /10 ) 
                    randomROIs<- sample(clean_ROI_list,ROIs_to_sample)
                    print(paste0("Total random ROIs selected is : ", length(randomROIs)) )
                    tag_var_list = c('H','F')
                    #### do dual_record stuff
#                     for(i in 1:total_ROIs) {
#                             current_tag = "H"#tag_var_list[i]
#                             #subset_traces = traces_df_fixed %>% dplyr::filter(trackROI_key %in% randomROIs[i])
#                             #tmp_plot<-plot_dualRecords(subsetted_traces = subset_traces,interspike_threshold = interSpike_thresh,tag = current_tag)
#                             if(i == 1){
# #### save to global env
#                              #  traceE  <<- tmp_plot

#                             }
#                             if(i == 2){
# #### save to global env

#                               # traceC  <<- tmp_plot
#                             }
#                         }

                    transmission_peak_groupers =  c('chemical_condition','vid_key','trackROI_key','ROINumber','peakID')
                    transmission_groupers = c('chemical_condition','vid_key','trackROI_key','ROINumber')
                    
                    clean_transmissionCheck <- traces_df_fixed %>% 
                                      dplyr::filter(trackROI_key %in% clean_ROI_list) %>%
                                      dplyr::filter(peakID != "NotPeak") %>% 
                                      ungroup() %>%
                                      group_by_at(transmission_peak_groupers) %>% 
                                      slice(which.max(dFF)) %>% 
                                      ungroup() %>% 
                                      group_by_at(transmission_groupers) %>% 
                                      arrange(chemical_condition,trackROI_key,absoluteTime) %>%
                                      mutate(
                                             interSpike = absoluteTime - lag(absoluteTime, default = first(absoluteTime)),
                                             isValid = ifelse(sensor == "JF646" & lag(sensor) == "GluSnFR3", TRUE, FALSE),
                                             timeCheck = ifelse(isValid == TRUE & interSpike < interSpike_thresh, TRUE, FALSE),
                                              isTransmission =  case_when(sensor=="GluSnFR3" & lead(sensor) == "JF646" ~lead(timeCheck), 
                                                                          sensor == "GluSnFR3" & is.na(lead(sensor)) ~ FALSE,
                                                                          sensor == "GluSnFR3" & lead(sensor) == "GluSnFR3" ~ FALSE),
                                              interSpike_fixed = ifelse(isTransmission == TRUE, lead(interSpike), NA),
                                              checkLonely = ifelse(sensor=="GluSnFR3" & is.na(lead(sensor)), TRUE, FALSE),
                                              Gluamplitude = ifelse(sensor == "GluSnFR3" & isTransmission == TRUE, dFF, NA),
                                               JFamplitude = ifelse(sensor == "GluSnFR3" & isTransmission == TRUE, lead(dFF), NA)) 




                    get_clean_transmission_ROIs <- clean_transmissionCheck %>% dplyr::filter(!is.na(isTransmission) , 
                                                                                                chemical_condition == "ctrl")
                    clean_transmission_ROIs <- unique(get_clean_transmission_ROIs$trackROI_key)


                    clean_transmission_analysis<- clean_transmissionCheck %>% dplyr::filter(trackROI_key %in% clean_transmission_ROIs) %>% ungroup() %>% group_by(chemical_condition, sensor, isTransmission) %>% tally()
                    clean_transmission_summary<- clean_transmission_analysis %>% mutate(isTransmission_fixed = case_when(isTransmission == TRUE ~ "success",
                                                                                                                            isTransmission == FALSE ~ "failure",
                                                                                                                            sensor == "JF646" & is.na(isTransmission)~ "orphan"
                                                                                                                            )

                                                                                        ) %>%
                                                                                group_by(chemical_condition) %>%
                                                                                summarise(#n_ROIs = n_distinct(trackROI_key),
                                                                                            n_success = n[which(isTransmission_fixed=="success")],
                                                                                            n_failure = n[which(isTransmission_fixed=="failure")],
                                                                                            n_orphan = n[which(isTransmission_fixed=="orphan")],
                                                                                            n_attempt = n_success+n_failure,
                                                                                            transmission_rate = n_success/n_attempt,
                                                                                            orphan_rate = n_orphan/(n_orphan+n_attempt))



                    clean_transmissionCheck_forplot  <- clean_transmissionCheck %>% group_by_at(transmission_groupers) %>%
                                       dplyr::filter(sensor == "GluSnFR3",chemical_condition == "ctrl") %>%
                                       mutate(deltaT  = absoluteTime - lag(absoluteTime, default = NA),
                                                Gluamplitude = dFF ,
                                                isTransmission_fixed = case_when(isTransmission == TRUE ~ "Transmission Success",
                                                                                    isTransmission == FALSE ~ "Transmission Failure")) %>%
                                       ungroup() %>%
                                       group_by(isTransmission_fixed) %>%
                                       dplyr::filter(isTransmission_fixed %in% c("Transmission Success","Transmission Failure")) %>%
                                       mutate(med_Gluamplitude = median(as.numeric(Gluamplitude)),
                                                med_JFamplitude = median(as.numeric(JFamplitude)),
                                                )





                  transmission_rateCheck <- clean_transmissionCheck  %>%
                                       group_by(chemical_condition,vid_key,trackROI_key,ROINumber) %>%
                                       dplyr::filter(sensor == "GluSnFR3") %>%
                                       mutate(deltaT  = absoluteTime - lag(absoluteTime, default = NA) ) %>% ### double check this implementation
                                       summarise(amplitudeJF_success = mean(JFamplitude, na.rm=TRUE),
                                                 med_Glu_amplitude = median(dFF,na.rm=TRUE),
                                                 med_JF_amplitude = median(JFamplitude,na.rm=TRUE),
                                                 amplitudeGlu_failure = mean(dFF[which(isTransmission == FALSE)],na.rm=TRUE),
                                                 amplitudeGlu_success = mean(dFF[which(isTransmission == TRUE)],na.rm=TRUE),
                                                 deltaT_failure = mean(deltaT[which(isTransmission == FALSE)],na.rm=TRUE),
                                                 deltaT_success = mean(deltaT[which(isTransmission == TRUE)],na.rm=TRUE),
                                                 num_JF = n_distinct(JFamplitude),
                                                 JF_freq = num_JF / 22.91, #seconds
                                                 num_transmission_fail = n_distinct(peakID[which(isTransmission == FALSE)]),
                                                 num_transmission_success = n_distinct(peakID[which(isTransmission == TRUE)]),
                                                  interSpikeCheck = mean(interSpike_fixed[which(isTransmission == TRUE)],na.rm=TRUE),
                                                  num_transmission_attempt = num_transmission_success+num_transmission_fail,
                                                  transmission_rate = num_transmission_success/num_transmission_attempt  ) %>%
                                       mutate(exchange_status = case_when(chemical_condition == "ctrl" ~ "Before",
                                                                                            chemical_condition == "vehicle" | chemical_condition == "APV" ~ "After") )%>%
                                                         
                                       ungroup() 
                    #remove_zeros_after<- transmission_rateCheck %>% dplyr::filter(exchange_status == "After")#,num_transmission_attempt >= 1)#, num_transmission_success >=0)
                    remove_zeros_before<- transmission_rateCheck %>% dplyr::filter(exchange_status == "Before",num_transmission_success >=1)
                    #variance_filter<-
                    #transmission_rateCheck <- transmission_rateCheck %>% dplyr::filter(trackROI_key %in% unique(remove_zeros_after$trackROI_key))
                    transmission_rateCheck <- transmission_rateCheck %>% dplyr::filter(trackROI_key %in% unique(remove_zeros_before$trackROI_key))


#sample ROIs

                    clean_transmission_ROIs<- transmission_rateCheck %>% dplyr::filter(trackROI_key %in% ROIs_to_save)
                    clean_transmission_ROIs_list = unique(clean_transmission_ROIs$trackROI_key)
                    total_ROIs<- length(clean_transmission_ROIs_list)
                    #ROIs_to_sample = ceiling( total_ROIs /10 ) 
                    #randomROIs<- sample(clean_ROI_list,ROIs_to_sample)
                    print(paste0("Total random ROIs selected is : ", length(clean_transmission_ROIs_list)) )
                    tag_var_list = c('I','G')
                    #### do dual_record stuff
#                     for(i in 1:total_ROIs) {
#                             current_tag = tag_var_list[i]
#                             subset_traces = traces_df_fixed %>% dplyr::filter(trackROI_key %in% clean_transmission_ROIs_list[i])
#                             tmp_plot<-plot_dualRecords(subsetted_traces = subset_traces,interspike_threshold = interSpike_thresh,tag = current_tag)
#                             if(i == 1){
# #### save to global env
#                              traceI  <<- tmp_plot

#                             }
#                             if(i == 2){
# #### save to global env

#                               traceG  <<- tmp_plot
#                             }
#                         }
                     
levels = c("ctrl","vehicle","APV")
labels = c("Pre-Treatment","Vehicle","AP5")
                    transmission_rateCheck$chemical_condition <- factor(transmission_rateCheck$chemical_condition, labels=labels,levels= levels)
                    transmission_rateCheck$binary_chemical<-ifelse(transmission_rateCheck$chemical_condition=="Pre-Treatment",1,2)
                    transmission_rateCheck$groupfact<- factor(transmission_rateCheck$binary_chemical)
                    transmission_label_vehicle <- clean_transmissionCheck %>% dplyr::filter(chemical_condition == "APV") %>% ungroup() %>% mutate(groupID = "AP5") %>% select(trackROI_key,groupID)
                    transmission_label_APV <- clean_transmissionCheck %>% dplyr::filter(chemical_condition == "vehicle") %>% ungroup() %>% mutate(groupID = "Vehicle") %>% select(trackROI_key,groupID) 
                    transmission_labels <- bind_rows(transmission_label_APV,transmission_label_vehicle)

                    transmission_rate_fixed<- left_join(transmission_rateCheck, transmission_labels) %>% dplyr::filter(!is.na(groupID) )

                    levels_groupID = c('Vehicle','AP5')
                    transmission_rate_fixed$groupID=factor(transmission_rate_fixed$groupID, levels= levels_groupID)



#             ggsave(filename=paste0("basePlot_transmissionPlot.png"),plot=transmissionPlot, device="png",dpi=600, units="in",width=lineplot_width,height=plot_height)



                    
                    count_veh_transmission = transmission_rate_fixed %>% dplyr::filter(chemical_condition == "Vehicle")
                    print(paste0("Number of tracked ROIs for vehicle condition synaptic transmission is : ", length(unique(count_veh_transmission$trackROI_key))) )

                    count_APV_transmission = transmission_rate_fixed %>% dplyr::filter(chemical_condition == "AP5")
                    print(paste0("Number of tracked ROIs for AP5 condition synaptic transmission is : ", length(unique(count_APV_transmission$trackROI_key))) )
                    

                    labels_revised = c('Before', 'After')

                    tmp_transmission_rate_fixed <- transmission_rate_fixed #%>% dplyr::filter(groupID == "AP5")
                    transmission_rate_v1 <- ggplot(tmp_transmission_rate_fixed, aes(x = binary_chemical, y = num_transmission_fail, alpha=groupID) ) +
                                            geom_line(aes(group =trackROI_key), alpha=0.1, colour="grey21", linewidth=1)+
                                           geom_point(aes(group = trackROI_key), alpha=0.1, colour="grey21", size=2)+
                       
                                           stat_summary(aes(group=groupID, colour=groupID,linetype=groupID), geom="line", fun=mean,linewidth=2,alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID), geom="errorbar", fun.data=mean_se,width=0.1,size=2, alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID,shape=groupID), geom="point", stroke=3,fill='white',fun=mean,size=5,alpha=1,na.rm=TRUE)+
                                           labs( x="",
                                                  y=expression(N[transmission_fail]),#expression(P[transmission]),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="")+
                                           scale_colour_manual(name = "", values = c("Vehicle" = "black", 'AP5' = 'red'))+
                                           scale_shape_manual(name = "", values = c("Vehicle" = 16, 'AP5' = 21))+
                                           scale_linetype_manual(name = "", values = c("Vehicle" = "solid", 'AP5' = 'solid'))+
                                           scale_alpha_manual(name="",values = c("Vehicle" = 1.0, "AP5" = 1.0))+
                                           
                                           # coord_cartesian(xlim=c(0.9,2.1),ylim=c(0.5,3.5))+
                                           #scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
                                           scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~groupID)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(#axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank()
                                                 #legend.position=c(0.5,1.05),
                                                 #legend.justification="center"
                                            )

                                            #ggsave(filename=paste0("transmission_rate_v1.png"),plot=transmission_rate_v1, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                            save_plot(filename=paste0("transmission_rate__as_num_failures_v1"),plot=transmission_rate_v1, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)
                    
                    levels(tmp_transmission_rate_fixed$groupID) = c("Vehicle","AP5")
                    tmp_transmission_rate_fixed$groupID = factor(tmp_transmission_rate_fixed$groupID, levels=c("Vehicle","AP5"),labels=c("Control", expression("100"~mu*"M AP5")))

                     transmission_rate_v1 <- ggplot(tmp_transmission_rate_fixed, aes(x = binary_chemical, y = num_transmission_success, alpha=groupID) ) +
                                            geom_line(aes(group =trackROI_key), alpha=0.1, colour="grey21", linewidth=1)+
                                           geom_point(aes(group = trackROI_key), alpha=0.1, colour="grey21", size=2)+
                       
                                           stat_summary(aes(group=groupID, colour=groupID,linetype=groupID), geom="line", fun=mean,linewidth=2,alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID), geom="errorbar", fun.data=mean_se,width=0.1,size=2, alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID,shape=groupID), geom="point", stroke=2,fill='white',fun=mean,size=4,alpha=1,na.rm=TRUE)+
                                           labs( x="",
                                                  y=expression(N[transmission_success]),#expression(P[transmission]),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="I")+
                                           scale_colour_manual(name = "", values = c("black", 'red'))+#c("Vehicle" = "black", 'AP5' = 'red'))+
                                           scale_shape_manual(name = "", values = c(16, 21))+#c("Vehicle" = 16, 'AP5' = 21))+
                                           scale_linetype_manual(name = "", values = c("solid", 'solid'))+#c("Vehicle" = "solid", 'AP5' = 'solid'))+
                                           scale_alpha_manual(name="",values = c(1.0, 1.0))+#c("Vehicle" = 1.0, "AP5" = 1.0))+
                                           
                                           # coord_cartesian(xlim=c(0.9,2.1),ylim=c(0.5,3.5))+
                                           #scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
                                           scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~groupID,labeller=label_parsed)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_text(colour="black", size=scalefactor*28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_text(colour="black", size=24, family="sans"),
                                                 legend.position="none"
                                                 #legend.justification="center"
                                            )

                                            #ggsave(filename=paste0("transmission_rate_v1.png"),plot=transmission_rate_v1, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                            save_plot(filename=paste0("transmission_rate__as_num_success_v1"),plot=transmission_rate_v1, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)
                    peak_transmission_I <<- transmission_rate_v1

                     transmission_rate_v1 <- ggplot(tmp_transmission_rate_fixed, aes(x = binary_chemical, y = num_transmission_attempt, alpha=groupID) ) +
                                            geom_line(aes(group =trackROI_key), alpha=0.1, colour="grey21", linewidth=1)+
                                           geom_point(aes(group = trackROI_key), alpha=0.1, colour="grey21", size=2)+
                       
                                           stat_summary(aes(group=groupID, colour=groupID,linetype=groupID), geom="line", fun=mean,linewidth=2,alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID), geom="errorbar", fun.data=mean_se,width=0.1,size=2, alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID,shape=groupID), geom="point", stroke=3,fill='white',fun=mean,size=5,alpha=1,na.rm=TRUE)+
                                           labs( x="",
                                                  y=expression(N[transmission_attempt]),#expression(P[transmission]),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="")+
                                           scale_colour_manual(name = "", values = c("Vehicle" = "black", 'AP5' = 'red'))+
                                           scale_shape_manual(name = "", values = c("Vehicle" = 16, 'AP5' = 21))+
                                           scale_linetype_manual(name = "", values = c("Vehicle" = "solid", 'AP5' = 'solid'))+
                                           scale_alpha_manual(name="",values = c("Vehicle" = 1.0, "AP5" = 1.0))+
                                           
                                           # coord_cartesian(xlim=c(0.9,2.1),ylim=c(0.5,3.5))+
                                           #scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
                                           scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~groupID)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(#axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank()
                                                 #legend.position=c(0.5,1.05),
                                                 #legend.justification="center"
                                            )

                                            #ggsave(filename=paste0("transmission_rate_v1.png"),plot=transmission_rate_v1, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                            save_plot(filename=paste0("transmission_rate__as_num_attempt_v1"),plot=transmission_rate_v1, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)


                     transmission_rate_v1 <- ggplot(tmp_transmission_rate_fixed, aes(x = binary_chemical, y = transmission_rate, alpha=groupID) ) +
                                            geom_line(aes(group =trackROI_key), alpha=0.1, colour="grey21", linewidth=1)+
                                           geom_point(aes(group = trackROI_key), alpha=0.1, colour="grey21", size=2)+
                       
                                           stat_summary(aes(group=groupID, colour=groupID,linetype=groupID), geom="line", fun=mean,linewidth=2,alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID), geom="errorbar", fun.data=mean_se,width=0.1,size=2, alpha=1,na.rm=TRUE)+
                                           stat_summary(aes(group=groupID, colour=groupID,shape=groupID), geom="point", stroke=3,fill='white',fun=mean,size=5,alpha=1,na.rm=TRUE)+
                                           labs( x="",
                                                  y=expression(P[transmission_success]),#expression(P[transmission]),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="")+
                                           scale_colour_manual(name = "", values = c("Vehicle" = "black", 'AP5' = 'red'))+
                                           scale_shape_manual(name = "", values = c("Vehicle" = 16, 'AP5' = 21))+
                                           scale_linetype_manual(name = "", values = c("Vehicle" = "solid", 'AP5' = 'solid'))+
                                           scale_alpha_manual(name="",values = c("Vehicle" = 1.0, "AP5" = 1.0))+
                                           
                                           # coord_cartesian(xlim=c(0.9,2.1),ylim=c(0.5,3.5))+
                                           #scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))+
                                           scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           facet_grid(~groupID)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(#axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 strip.text.x=element_blank()
                                                 #legend.position=c(0.5,1.05),
                                                 #legend.justification="center"
                                            )

                                            #ggsave(filename=paste0("transmission_rate_v1.png"),plot=transmission_rate_v1, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                            save_plot(filename=paste0("transmission_rate__as_ratio_transmission_v1"),plot=transmission_rate_v1, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)




                    #labels_revised = c('Before', 'After')
                    check_interSpikes<- clean_transmissionCheck_forplot %>% dplyr::filter(isTransmission_fixed == "Transmission Success")#, interSpike_fixed>0)
                    transmission_interSpike_v1 <- ggplot(check_interSpikes, aes(x = isTransmission_fixed, y = interSpike_fixed*1000, colour=isTransmission_fixed) ) +
                                           geom_sina(size=4,alpha=0.7)+
                                           geom_violin(fill=NA, size=1,colour="black",draw_quantiles = c(0.25, 0.5, 0.75))+
                                           labs( x="",
                                                  y=expression(Delta*"t between paired events (ms)"),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="")+
                                           scale_colour_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           #scale_fill_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           #scale_alpha_manual(name="",values = c("Transmission Success" = 0.8,"Transmission Failure" = 0.5))+
                                           
                                           # coord_cartesian(xlim=c(0.9,2.1),ylim=c(0.5,3.5))+
                                           scale_y_continuous(limits=c(0,400))+
                                           #scale_x_continuous(breaks = c(1,2), labels = labels_revised)+#scale_y_continuous(limits=c())
                                           #facet_grid(~groupID)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(axis.text.x=element_blank(),#element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 axis.line.x = element_blank(),
                                                 axis.ticks.x= element_blank(),
                                                 strip.text.x=element_blank(),
                                                 legend.position = "none"
                                            )

                                            #ggsave(filename=paste0("transmission_interspike_v1.png"),plot=transmission_interSpike_v1, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                            save_plot(filename=paste0("transmission_interspike_v1"),plot=transmission_interSpike_v1, plot_width=plot_dim*0.5,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

                                            transmission_interspike_plot_E  <<- transmission_interSpike_v1




#### save to global env
                              
                    #transmission_rate_plot_K  <<- transmission_rate_v1


                    clean_transmissionCheck_scatterplot_data<- clean_transmissionCheck_forplot %>% dplyr::filter(isTransmission_fixed == "Transmission Success")

            amplitude_scatterplot<-ggplot(clean_transmissionCheck_forplot, aes(x=Gluamplitude, y=JFamplitude, group=isTransmission_fixed,colour=isTransmission_fixed))+
                                                geom_point(size=4,alpha=0.9)+
                                                geom_smooth(method = "lm",formula = y~ x, se = FALSE, colour = "black", size=2)+
                                                stat_cor( color = "black", geom = "text",label.x = 1.5,label.y=0.95,size=8)+ #aes(label = after_stat(rr.label)),
                                                stat_regline_equation(label.x = 1.5, label.y = 0.85,color="black",size=8) +

                                                
                                      labs( x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
                                                  y=expression(JF[646]*" amplitude ("*Delta*"F/F)"),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="C")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           coord_cartesian(xlim=c(0,3),ylim=c(0,1))+
                                           scale_y_continuous(breaks=c(0,0.5,1.0))+
                                           scale_x_continuous(limits=c(0,3),breaks=c(0,1,2,3))+#scale_y_continuous(limits=c())
                                           #facet_grid(isTransmission_fixed~.)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(legend.position = "none"
                                                 #axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 #strip.text.x=element_blank(),
                                                 #strip.text.y=element_blank(),#element_text(angle = 0, size=28, family="sans"),
                                                 #panel.spacing = unit(1, "cm"),
                                                 #legend.position = c(.95, .95),
                                                 #legend.justification = c("right", "top")
                                            )
                                 #ggsave(filename=paste0("failure_vs_success_scatterplot.png"),plot=amplitude_scatterplot, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                 save_plot(filename=paste0("failure_vs_success_scatterplot_v1"),plot=amplitude_scatterplot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

            scatter_amplitude_corr_C<<-amplitude_scatterplot





            amplitude_scatterplot<-ggplot(clean_transmissionCheck_forplot, aes(x=Gluamplitude, y=interSpike_fixed*1000, group=isTransmission_fixed,colour=isTransmission_fixed))+
                                                geom_point(size=4,alpha=0.9)+
                                                geom_smooth(method = "lm", formula=y~x,se = FALSE, colour = "black", size=2)+
                                                 stat_cor(color = "black", geom = "text",label.x = 1.5,label.y=0.95*500,size=8)+ #aes(label = after_stat(rr.label)), 
                                                stat_regline_equation(label.x = 1.5, label.y = 0.85*500,color="black",size=8) +

                                                
                                      labs( x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
                                                  y=expression(Delta*"t to paired "*JF[646]~"event (ms)"),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="D")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           #coord_cartesian(xlim=c(0,3),ylim=c(0,1))+
                                           #scale_y_continuous(breaks=c(0,0.5,1.0))+
                                           scale_x_continuous(limits=c(0,3),breaks=c(0,1,2,3))+#scale_y_continuous(limits=c())
                                           #facet_grid(isTransmission_fixed~.)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(legend.position = "none"
                                                 #axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 #strip.text.x=element_blank(),
                                                 #strip.text.y=element_blank(),#element_text(angle = 0, size=28, family="sans"),
                                                 #panel.spacing = unit(1, "cm"),
                                                 #legend.position = c(.95, .95),
                                                 #legend.justification = c("right", "top")
                                            )
                                 #ggsave(filename=paste0("failure_vs_success_scatterplot.png"),plot=amplitude_scatterplot, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                 save_plot(filename=paste0("dT_vs_Gluamp_scatterplot_v1"),plot=amplitude_scatterplot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

            scatter_Gluamplitude_dt_D <<-amplitude_scatterplot
             amplitude_scatterplot<-ggplot(clean_transmissionCheck_forplot, aes(x=interSpike_fixed*1000, y=JFamplitude, group=isTransmission_fixed,colour=isTransmission_fixed))+
                                                geom_point(size=4,alpha=0.9)+
                                                geom_smooth(method = "lm", se = FALSE, colour = "black", size=2)+
                                                
                                      labs( y=expression(JF[646]*" amplitude ("*Delta*"F/F)"),
                                                  x=expression(Delta*"t from paired iGluSnFR3 event (ms)"),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("Transmission Success" = "#FF33FF","Transmission Failure" = "grey21"))+
                                           #coord_cartesian(xlim=c(0,3),ylim=c(0,1))+
                                           scale_y_continuous(breaks=c(0,0.5,1.0))+
                                           #scale_x_continuous(limits=c(0,3),breaks=c(0,1,2,3))+#scale_y_continuous(limits=c())
                                           #facet_grid(isTransmission_fixed~.)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(legend.position = "none"
                                                 #axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 #strip.text.x=element_blank(),
                                                 #strip.text.y=element_blank(),#element_text(angle = 0, size=28, family="sans"),
                                                 #panel.spacing = unit(1, "cm"),
                                                 #legend.position = c(.95, .95),
                                                 #legend.justification = c("right", "top")
                                            )
                                 #ggsave(filename=paste0("failure_vs_success_scatterplot.png"),plot=amplitude_scatterplot, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                 save_plot(filename=paste0("JFamp_vs_dT_scatterplot_v1"),plot=amplitude_scatterplot, plot_width=plot_dim,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)

####STATISTICS
                                 A <- clean_transmissionCheck_forplot %>% dplyr::filter(isTransmission_fixed == "Transmission Failure")
                                 B <- clean_transmissionCheck_forplot %>% dplyr::filter(isTransmission_fixed == "Transmission Success")

                                 tmp_A <- A %>% ungroup() %>% select(Gluamplitude)
                                 tmp_A <- as.matrix(tmp_A)
                                 tmp_B <- B %>% ungroup() %>% select(Gluamplitude)
                                 tmp_B <- as.matrix(tmp_B)
                                 
    
                                run_stats <- function(parent_A,comparison_A, parent_B,comparison_B){#,list_store) {
                                    print(paste0("wilcox test"))
                                    test<-wilcox.test(comparison_A,comparison_B,paired=FALSE)
                                    print(test)
                                    comparison = paste(unique(parent_A$isTransmission_fixed),unique(parent_B$isTransmission_fixed),sep='-')
                                    data_wilcox<- data.frame(comparison = comparison,test = 'wilcoxon',p_value = test$p.value)
                                    print(paste0("ks.test for"))
                                    test<-ks.test(comparison_A,comparison_B)
                                    print(test)
                                    data_ks<- data.frame(comparison = comparison,test = 'ks_test',p_value = test$p.value)
                                    data_bind <- bind_rows(data_wilcox,data_ks)

                                    data_bind
                                }

                                print(paste0("Generating a statistical comparison for Gluamplitude between Transmission Success and Failure ") )
                                check_output<-run_stats(A,tmp_A,B,tmp_B)#,p_list)
    



            amplitude_histogram<-ggplot(clean_transmissionCheck_forplot, aes(x=Gluamplitude, group=isTransmission_fixed,colour=isTransmission_fixed,fill=isTransmission_fixed,alpha=isTransmission_fixed))+
                                                #geom_histogram(data=clean_noise, aes(x=noise, y=after_stat(ncount),group=groupID),fill='grey21',binwidth = binwidth, alpha=0.5,boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                #geom_freqpoly(data=clean_noise, aes(x=noise,y=after_stat(ncount)),colour='grey21',binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5)+
                                                geom_freqpoly(aes(y=after_stat(ncount)),binwidth = binwidth, alpha=0.8,,boundary =0,linewidth=1.5,show_guide = FALSE)+
                                                
                                                geom_histogram(aes(y=after_stat(ncount)),binwidth = binwidth, boundary = 0,position="identity",colour="black",linewidth=0.8)+
                                                
                                                geom_vline(aes(xintercept = med_Gluamplitude), colour= "black", lty = "dashed", linewidth = 2,alpha=1)+
                                                #stat_compare_means()+
                                      labs( x=expression("iGluSnFR3 amplitude ("*Delta*"F/F)"),
                                                  y=expression(N/N[max]),
                                                  ##title="",
                                                  #subtitle="Averaged response from three trials per protocol",
                                                  tag="E")+                                                      #"GluSnFR3-Vehicle" = "#336600", 'GluSnFR3-AP5' = '#339900', "JF646-Vehicle" ="#990066","JF646-AP5" ="#FF33FF"
                                           scale_colour_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           scale_fill_manual(name = "", values = c("Transmission Success" = '#339900',"Transmission Failure" = "grey21"))+
                                           scale_alpha_manual(name="",values = c("Transmission Success" = 0.8,"Transmission Failure" = 0.5))+
                                            coord_cartesian(xlim=c(0,3),ylim=c(0,1.1))+
                                           scale_y_continuous(expand=c(0,0), breaks=c(0,0.5,1.0))+
                                           scale_x_continuous(limits=c(0,3))+#scale_y_continuous(limits=c())
                                           facet_grid(isTransmission_fixed~.)+
                                           theme_tufte()+
                                           my.theme +
                                           theme(#axis.text.x=element_text(colour="black", size=28, family="sans", angle=45, hjust=1),
                                                 #strip.text.x=element_blank(),
                                                 strip.text.y=element_text(colour="black", size=24, family="sans",angle=0),
                                                 panel.spacing = unit(1, "cm"),
                                                 legend.position = "none",
                                                 legend.justification = c("right", "top")
                                            )
                                 #ggsave(filename=paste0("failure_vs_scuess.png"),plot=amplitude_histogram, device="png",dpi=600, bg="white",units="in",width=lineplot_width,height=plot_height)
                                 save_plot(filename=paste0("failure_vs_success_histogram_v1"),plot=amplitude_histogram, plot_width=plot_dim*1.5,plot_height=plot_dim, scale_f = scalefactor,dpi_val = 600)


#### save to global env
                                #scatter_success_C  <<- amplitude_scatterplot
                                histo_success_E  <<- amplitude_histogram









# ###### GET SPATIAL DATA ##### ROI OVERLAY + ARROW

# # ###### for loop for plotting spatial data
#         png_dirs = c(#"Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_3",
#                         "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_6",
#                         "Y:\\Sam/paper1_datasets/synapticTransmission_v1/APV/Figures/synaptic_transmission_v4/synTransmission_v1/forMapping/_registered_11")
#         prefix_str = "AVG_GluSnFR3_JF646_4Ca_cis_"
#         suffix_str = "_ctrl_repl01__FLATTENED_COLOR.png"
#         all_lut_str = "_ctrl_repl01__FLATTENED_COLOR.png"
#         dir_regex = "registered_(.*?)_"
#         save_bool=TRUE
#         #print(paste0("Available vid_keys to save out are: ", unique(traces_df$vid_key) ) ) 
#         #subset_df<- traces_df_fixed %>% dplyr::filter(vid_key %in% vid_keys_to_save, trackROI_key %in% clean_ROI_list)
#         #print(paste0("Column names are: ", colnames(subset_df)) ) 
#         #print(paste0("Vid_key is: ", unique(subset_df$vid_key) ) )
#         #print(unique(subset_df$protocol))
#         #vid_key_vec <- unique(subset_df$vid_key)
        
#         print("Taking a crack at generating png plots")
#         subset_df<- transmission_rateCheck %>% dplyr::filter(vid_key %in% vid_keys_to_save,chemical_condition == "Pre-Treatment")

#         #ROInames <- str_extract(ROIs_to_save,"ROI\\d\\d\\d\\d") 
#         vid_key_vec <- unique(subset_df$vid_key)
#         vars_to_plot<- c('med_Glu_amplitude', "med_JF_amplitude",'JF_freq','transmission_rate')#c("releaseProbability_perROI")
#         save_bool=TRUE
#         guide_bool = TRUE

#         guide_limits = list(c(0.5,2),c(0.1,0.8),c(0,0.25),c(0,1))
#         guide_titles = c(expression(iGlu[Delta*F/F]),expression(JF[646][Delta*F/F]),expression(JF[646]~(Hz)),expression(P[transmission]))# # expression(CV[tau*"decay"])#
        
#         #expression("Avg."~Evoked[Delta*"F/F"] / Spont[Delta*"F/F"])
#         tmp_tag_var_override = c("F","G", "H")


#         ###scale bar params! set the left-corner 
#         user_x = 51.2-10.6
#         user_y = 50.2
#         tag_vars = c("F","H"    )



#         tag_vars = c("H","F")
#         print("getting ROInames") 
#         ROInames <- str_extract(save_ROIs,"ROI\\d\\d\\d\\d")   
#         print(ROInames)
#         for(j in 1:length(png_dirs)){
#             current_png_dir = png_dirs[j]
#             print(paste0("current png directory is: ",current_png_dir) )
#             print(paste0("Vid_key is : ",vid_key_vec[j]) )
               
#             tmp_subset_transmission <- subset_df %>% dplyr::filter(vid_key == vid_key_vec[j]) %>% mutate(fix_Ca = "4Ca")
#             print(head(tmp_subset_transmission))

#             for(k in 1:length(vars_to_plot)){            
#                                 guide_bool = TRUE
#                                 var_to_plot = vars_to_plot[k]
#                                 tmp_guide_title = guide_titles[k]
#                                 tmp_limits = guide_limits[[k]]
#                                 tmp_tag = ""#tmp_tag_var_override[j]
#                                 if(k == 1){tmp_scalebar = TRUE} else {tmp_scalebar = FALSE}
#                                 # tmp_plot = png_plotter(df=tmp_subset_transmission,
#                                 #                                     png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,
#                                 #                                     vars_to_plot=var_to_plot, 
#                                 #                                     save_bool=save_bool,guide_bool=guide_bool,tag_var =tmp_tag,
#                                 #                                     ROI_IDs = ROInames,guide_title=tmp_guide_title,guide_lim=tmp_limits,user_x = user_x, user_y = user_y,scalebar_switch=tmp_scalebar,add_arrows=FALSE,
#                                 #                                     scalebar_length=10,binfactor=2)

#             }
#         }

#         for (k in 1:total_ROIs){
#                 subset_traces = traces_df_fixed %>% dplyr::filter(trackROI_key %in% clean_transmission_ROIs_list[k])
#                 current_ROIs = unique(subset_traces$ROINumber)            
                
#                current_png_dir = png_dirs[k]
#                current_tag = tag_vars[k]
#                save_bool=TRUE    
                    
#                        png_plot<-png_plotter_ROIoverlay(df=subset_traces,png_dir=current_png_dir,prefix_str=prefix_str,suffix_str=suffix_str,dir_regex=dir_regex,save_bool=save_bool,tag_var = current_tag,ROI_IDs = current_ROIs)
                        

#                         if(k == 1){
# #### save to global env
#                             regionH  <<- png_plot


#                             }
#                         if(i == 2){
# #### save to global env
#                            regionF  <<- png_plot

#                             }         
                         
                    

#         }


                        





list.output = list(clean_ROI_list,
                    get_peak_maxima,
                    clean_transmissionCheck,
                    clean_transmission_analysis,
                    clean_transmission_summary,
                    count_peaks,
                    clean_transmission_ROIs,
                    transmission_rateCheck,
                    transmission_rate_fixed,
                    transmission_labels)#,
                    #remove_zeros)
list.output
}





