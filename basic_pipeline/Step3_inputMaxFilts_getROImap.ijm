// This script recurses into every folder in a subdirectory and looks for a user-defined suffix
//When it finds the suffix, it runs a Find Maxima... 
//Analyze particles...
//Enlarge
//and saves out the ROIs as a separate zip file (here, ROIs.zip) in a subdirectory named after the input file + _ROImaps
//This macro intends to threshold only the top 1% of pixels. 
//User should identify reasonable threshold prior to running batch based on the dataset.




// @File(label = "Input directory", style = "directory") input
// @File(label = "Output directory", style = "directory") output
// @String(label = "File suffix", value = ".tif") suffix

/*
 * Macro template to process multiple images in a folder
 */

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
dx = 5;
//the percentage (e.g. 0.01 = 1%) above which to threshold the pixels of the maxIntensity Filtered Image. 

input = getDirectory("Please choose a primary directory to process"); //prompts user to input a base folder to analyze ("directory" above)
suffix = "maxFiltered.tiff" //change this line each run to conform to maxFiltered naming convention
extract_xy_coords = 1 // 1 if true, 0 if false
if (extract_xy_coords == 1 ) {
	run("Set Measurements...", "centroid center nan redirect=None decimal=3");
}

//suffix = ".tif";
setBatchMode(true);
processFolder_ROIs(input);
//processFolder_ROIs(input);
print("Finished getting the ROImaps");

beep();



// function to scan folders/subfolders/files to find files with correct suffix
function processFolder_ROIs(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder_ROIs("" + input + list[i]);
		if(endsWith(list[i], suffix))
			getROIs(input, list[i]);
	}
}

function getROIs(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + file);
	path = input + file;
	open(path);
	title=getTitle();
	//run("Brightness/Contrast...");
	//run("Enhance Contrast", "saturated=0.35");
	//run("Apply LUT");
			//getDimensions(width, height, channels, slices, frames);
			//totalpx = width * height;
			percentage = 0.05;
			getStatistics(area, mean, min, max, std, histogram);
			target = (1-percentage)*getWidth()*getHeight();
			nBins = max;
			getHistogram(values,counts,nBins);
			sum = 0; threshold = -1;
			for(i=0; i<nBins; i++){
    			sum += counts[i];
    			if( (sum >= target) & (threshold < 0) ){ threshold = i; }
				}
			setThreshold(threshold,nBins);
			setOption("BlackBackground", false);	
			run("Threshold...");	
           	run("Convert to Mask");
           	run("Despeckle");
			run("Dilate");
			//run("Dilate");
			//run("Dilate");
			tmpDir = File.getDirectory(path) + File.separator;			       	
			binarySave = d2s(percentage, 2) + "binarymask.tif";
           	
			//Get ROIs between 10-1000 pixels 
			run("Analyze Particles...", "size=100-10000 exclude add"); //Find list of ROIs with small size
	
			n = roiManager("count");
			if(n > 0) {
           	savefile = "ROImap.zip";
			output= tmpDir + savefile;
			roiManager("Save", output);
			print("Saving ROIs to:" + output);
				
				if(extract_xy_coords == 1) {
					coord_savefolder = "xy_coords" + File.separator;
					File.makeDirectory(tmpDir + coord_savefolder);
					roiManager("Measure");
					for (row = 0; row < nResults; row++) {
							setResult("Image Name", row, title);
						}
					updateResults();
					saveAs("Results",tmpDir  + coord_savefolder +  "xy-coords_Results.csv");
					}
			}
           	saveAs("Tiff", tmpDir + binarySave);

    		print("ROIs identified: " + n);
    		print("Threshold applied should have been: top "+percentage*100 + "% of all pixels, i.e. " + threshold + "," + nBins);   
	close(title);
	cleanup();
}



function cleanup() {
	
	//Close any image windows and results windows, reset the ROI manager
	list = getList("image.titles");
	if (list.length>0) {
		run("Close All");
	}
	if (isOpen("Results")) {
		close("Results");
	}
	roiManager("reset");
}

