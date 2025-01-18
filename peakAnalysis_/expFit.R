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

					expFit<- df %>%
							dplyr::filter(peakID != "NotPeak", absoluteTime < 1.55) %>% #unit test ###dataset from munc13 correlation, max time = 1.5
							group_by_at(groupers) %>%
							slice(which.max(dFF):length(dFF)) %>%
							mutate(count_obs = n()) %>%
							dplyr::filter(count_obs >= 3) %>%
							mutate(normTime = absoluteTime - min(absoluteTime))  %>%
							#dplyr::filter(normTime < 1.0) %>%
		    					ungroup() %>%
		    					group_by_at(groupers) %>%
		    					nest() %>%
		   					  	mutate(
		   						        fit = purrr::map(data, ~nls_multstart(dFF ~ A*exp(-normTime/tau_decay),
		   						        				data=.x,
		   						        				iter = 1000,
		   						        				start_lower = c(A = 0.15, tau_decay= 0.01 ),
		   						        				start_upper = c(A = 6, tau_decay= 0.3),
		   						        				supp_errors = 'Y',
		   						        				na.action = na.omit,
		   						        				lower = c(A = 0, tau_decay= 0))),
		   						       	tidied = map(fit, tidy),
					            		augment = map(fit, augment)

							)


					peekFit<-expFit%>%
			  			unnest(augment) 

		   			
			  		fitData<- peekFit
		        


	  	
		
		fitData	  
}
