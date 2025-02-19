# load packages just in case

#suppress warnings
oldw <- getOption("warn")
options(warn = -1)



library(devtools)
library(nls.multstart) 
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)


expDecayFit<- function(df, groupers) {			
			#tic()
					expFit<- df %>%
							dplyr::filter(peakID != "NotPeak") %>% #unit test
							group_by_at(groupers) %>%
							slice(which.max(dFF):length(dFF)) %>%
							dplyr::filter(length(dFF)>=2) %>%
							mutate(normTime = absoluteTime - min(absoluteTime))  %>%
		    					ungroup() %>%
		    					group_by_at(groupers) %>%
		    					nest() %>%
		   					  	mutate(
		   						        fit = purrr::map(data, ~nls_multstart(dFF ~ A*exp(-normTime/tau_decay),
		   						        				data=.x,
		   						        				iter = 1000,
		   						        				start_lower = c(A = 0, tau_decay= 0 ),
		   						        				start_upper = c(A = 10, tau_decay= 1),
		   						        				supp_errors = 'Y',
		   						        				na.action = na.omit,
		   						        				lower = c(A = 0, tau_decay= 0))),
		   						       	tidied = map(fit, tidy),
					            		augment = map(fit, augment)

							)


					peekFit<-expFit%>%
			  			unnest(augment) 

		   			
			  		fitData<- peekFit
		        


	  	#toc()
		#print("Finished fitting an exponential to the decay phase of each peak.")
		
		fitData	  
}
