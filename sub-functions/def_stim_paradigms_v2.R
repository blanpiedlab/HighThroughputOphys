#def_stim_paradigms_v2.R
#Pseudo-code - this one's gonna be tricky!
subDir <- "checkStimFinder"   


source("\\\\blanpiedserver/NASShare3/Sam/Sam scripts/Active scripts/pipeline v3/setFigureDirectory.R") 






def_stim_paradigms<- function(dataframe, keys) {


## Step 1. Generate keys for each subset of the dataframe. 
.keyvars = rlang::syms(keys)
stim = fitsLabeled %>% 
			dplyr::filter(str_detect(stimParadigm, 'Hz')) %>%
			mutate(stimKey = paste(!!!.keyvars, sep="-"),
					stimOnly = stimParadigm)
unique.stimKeys = unique(stim$stimKey)
stimList = list() #instantiate stimList


## Step 2. Generate values for each subset of the dataframe.
					
for (j in 1:length(unique.stimKeys)) {

			keyName = unique.stimKeys[[j]]
			
			j_video = stim %>% ungroup() %>% dplyr::filter(stimKey %in% keyName)
			
			
			stim.string = unique(j_video$stimParadigm)
			numeric.stim = as.numeric( gsub( "_Hz", "", as.matrix( stim.string ) ) ) 
			interStimulus = 1/numeric.stim
			
			
			interFrame = mean(j_video$interFrame, na.rm=TRUE)
			stim_delay = 0.5
			framePad = ceiling( ( stim_delay-min(j_video$absoluteTime) )/interFrame )
			lastFrame = ceiling(interStimulus/interFrame)+framePad
			
			time.min = min(j_video$absoluteTime)
			last.sync.time = 0.20+time.min
			
			
			
			
			
			
			peakMax = suppressMessages(j_video %>% ungroup() %>% group_by(stimKey,ROINumber) %>%
			                          dplyr::filter(absoluteTime < last.sync.time) %>%
			                          summarise(peakMax = max(dFF[which(peakID != "NotPeak")])) 
			                )
		  	maxROI = peakMax[which.max(peakMax$peakMax),]
		  		
  
  #first exit from the for loop : Did not find a finite peak maximum in the expected area of the trace.
			  	if (is.finite(maxROI$peakMax) == FALSE ) {
			    
					    print("No peak could be identified. Stimulus train set to default.")
					    

					    default.onset = 0.465
					    final.dur = default.onset+25
					    set.sequence = seq(from=default.onset, to = final.dur, by=interStimulus)
					    stimList[[keyName]] <- set.sequence

						print(paste0("Entered a stimList entry for element :", j))

								    
		    			next


  				} 
    
    			print("A peak maximum was found within the first stimulus window.")
    			if (maxROI$peakMax < 0.5 ) {
			    
					    print("The identified peak maximum did not exceed our threshold. Stimulus train set to default.")
					    

					    default.onset = 0.465
					    final.dur = default.onset+25
					    set.sequence = seq(from=default.onset, to = final.dur, by=interStimulus)
					    stimList[[keyName]] <- set.sequence

						print(paste0("Entered a stimList entry for element :", j))

								    
		    			next


  				} 
    			
    
			    k=3
			    j_video_max = j_video %>% 
			      dplyr::filter(ROINumber %in% maxROI$ROINumber) %>%
			      mutate(smoothed_dFF = rollapply(dFF,k,mean, align="center",fill=NA)) %>%
			      dplyr::filter(!row_number() %in% c(1, n()))
			                    
                    
                  # begin peak onset interpolation routine  
                    time = j_video_max$absoluteTime[1:lastFrame]
                    y = j_video_max$smoothed_dFF[1:lastFrame]
                    dFF.df = data.frame(time,y)
                    sync.range = dFF.df[which(dFF.df$time < last.sync.time),] 
                    max.dFF.time = sync.range$time[which.max(sync.range$y)]
                      
                      
                    
                    
                    dFF.spline = smooth.spline(time,y,df=200)
                    get.deriv = predict(dFF.spline,time,deriv=2)
                    
                    x.spline = dFF.spline$x
                    y.spline = dFF.spline$y
                    spline.df = data.frame(x.spline,y.spline)
                    
                    x.deriv = get.deriv$x
                    y.deriv = get.deriv$y
                    deriv.df = data.frame(x.deriv,y.deriv)
                    
                    
                    onset.range = deriv.df[which(deriv.df$x.deriv < max.dFF.time & deriv.df$x.deriv >= 0.44),]
                    
        #second exit from the for loop : the found position of the maximum.dFF is too close to our minimum onset time for interpolation, i.e. no x values exist after filtering.  
          
                if  (length(onset.range$x.deriv) ==0) {
                  
		                print("Found max.dFF is too close to minimum trace time. Stimulus train set to default.")
		                

		                default.onset = 0.465  
		                final.dur = default.onset+25 
		                set.sequence = seq(from=default.onset, to = final.dur, by=interStimulus)
		                stimList[[keyName]] <- set.sequence

						print(paste0("Entered a stimList entry for element :", j))

								    
						next

                  
                } 
                	
                	largest_rate_of_change = onset.range$x.deriv[ which( onset.range$y.deriv==max(onset.range$y.deriv) )  ]
                    find.onset = onset.range$x.deriv[which(onset.range$x.deriv==largest_rate_of_change)]
                    

        #third exit from the for loop : Our found onset time for the first stimulus is not null, is length > 0, but exceeds the upper bound of what we expected (0.475 usually overlaps with the rise phase of the first peak or the maximum)
                if  ( is.null(find.onset) | length(find.onset) == 0 | !is.null(find.onset) & length(find.onset) > 0 & find.onset >= 0.475 ) {
                  
		                print("Found onset seems suspect or broken. Stimulus train set to default.")
		                
		                default.onset = 0.465  
		                final.dur = default.onset+25
		                set.sequence = seq(from=default.onset, to = final.dur, by=interStimulus)
						stimList[[keyName]] <- set.sequence

						print(paste0("Entered a stimList entry for element :", j))

								    
						next

                  
     	#if we can get through all of these checks, generate plots which show us the behavior of the routine.             
	            	}



                    final.dur = find.onset+25
                    
                    set.sequence = seq(from=find.onset, to = final.dur, by=interStimulus)
                    
                    plotCheck<-ggplot(j_video, aes(x=absoluteTime, y=dFF))+
                      geom_path(aes(group=ROINumber), alpha=0.5,colour='grey',size=1)+
                      geom_line(data=j_video_max, aes(x=absoluteTime, y= smoothed_dFF), colour="red",size=2)+
                      geom_line(data=spline.df,aes(x=x.spline,y=y.spline), colour="blue", size=1)+
                      #geom_line(data=deriv.df,aes(x=x.deriv,y=y.deriv),colour="orange",size=1,lty="dashed")+
                      geom_vline(xintercept=c(set.sequence))+
                      coord_cartesian( xlim= c(time.min,final.dur/5) )+
                      labs(title=keyName)
                    
                    
                    
                    ggsave(filename=paste0("foundStimStart_",j,".png"),plot=plotCheck, device="png",dpi=300, units="in",width=12,height=12)
                    
                    
                    plotCheck<-ggplot(j_video, aes(x=absoluteTime, y=dFF))+
                      geom_path(aes(group=ROINumber), alpha=0.5,colour='grey',size=1)+
                      geom_line(data=j_video_max, aes(x=absoluteTime, y= smoothed_dFF), colour="red",size=2)+
                      geom_line(data=spline.df,aes(x=x.spline,y=y.spline), colour="blue", size=1)+
                      geom_point(data=deriv.df,aes(x=x.deriv,y=y.deriv),colour="orange",size=3)+
                      geom_vline(xintercept=c(set.sequence))+
                      coord_cartesian( xlim= c(time.min,final.dur/5) )+
                      labs(title=keyName)
                      
                 
                    
                    
                    ggsave(filename=paste0("foundStimStart_",j,"_show2ndDeriv.png"),plot=plotCheck, device="png",dpi=300, units="in",width=12,height=12)
                    
  					

			stimList[[keyName]] <- set.sequence
			print(paste0("Entered a stimList entry for element :", j))

			

}


			stimList



}

