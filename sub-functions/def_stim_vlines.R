#function to call in expFit_tracePlotter - defines which vector to choose for the value of stim_vector

def_stim_vlines<- function(x = tmp_stim_ID) {

		timeName = x # this is the same way we defined it in def_stim_paradigms.R

	if ( !is.null(timeName) &  (timeName %in% names(stimList) ) == TRUE) {
		
	
		stimList[[timeName]]
	
	} else {

		NULL # stim ID corresponds to sham Stimulus condition or something else undefined
	}


}
