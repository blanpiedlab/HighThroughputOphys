#function to call in expFit_tracePlotter - defines which vector to choose for the value of stim_vector

def_stim_vlines<- function(x = tmp_stim_ID) {

	if (x == "0.5 Hz") {
	
		halfHz
	
	} else if (x == "1 Hz") {
	
		oneHz
	
	} else if (x == "2 Hz") {
		
		twoHz
	
	} else {

		NULL # stim ID corresponds to sham Stimulus condition or something else undefined
	}


}
