## Flexibly define all available stimulus paradigms provided they conform to an expected string formula 

def_stim_paradigms<- function(dataframe) {
    
			stim = dataframe %>% dplyr::filter(str_detect(stimParadigm, 'Hz'))
			
			stimParadigms = unique(stim$stimParadigm)  
			stimList = list()
      
			for ( i in 1:length(stimParadigms) ) {
			
						
			      
			      timeName = paste0(stimParadigms[[i]],"_times")
						
						numeric.stimParadigms = as.numeric( gsub( "_Hz", "", as.matrix( stimParadigms[[i]] ) ) ) 
						
						stimTimes <- seq(from=0.49, to = 25.49, by = 1/numeric.stimParadigms)											
			      stimList[[timeName]] <- stimTimes

      }

  stimList
}

