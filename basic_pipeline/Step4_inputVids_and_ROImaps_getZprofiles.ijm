// This script recurses into every folder in a subdirectory and finds the ROImap.zip file. 
// Applies ROImap.zip to the video (suffix = "_.tif") in the folder and extracts each 
// z-profile as an intensity-time trace (.csv).
// z-profiles are stored in subfolders with the same name as the video they came from; 
// these subfolders are aggregated into //ROIoutputs, which is manipulated afterwards in R. 

dx = 5;

input = getDirectory("Please choose a primary directory to process"); //prompts user to input a base folder to analyze ("directory" above)

prefix = "C1-";
suffix = ".tif";
setBatchMode(false);
processFolder_zProfiles(input);
print("Finished getting the zProfiles.");

beep();



// function to scan folders/subfolders/files to find files with correct suffix
function processFolder_zProfiles(input) {
	list = getFileList(input);
	//count = 0;
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder_zProfiles("" + input + list[i]);
		if(startsWith(list[i], prefix) && endsWith(list[i], suffix))
			getZprofiles(input, list[i]);
			//count=count+1;
			
	}
	//print("CAN'T MISS THIS MESSAGE: COUNT OF FILES IS ", count);
}




function getZprofiles(input, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	//subdir = File.getDirectory(file);
	//if (File.exists(File.getNameWithoutExtension(input+subdir+file) + "maxFilteredROImap.zip") != 1) {print("     Skipping file " + filename + ": Missing ROImap.zip."); return;}
    //if(indexOf(file, "0.01binary.tif") != 1) {
    	
    
    path = input + file;
    ROIpath = File.getParent(path) + File.separator;
	ROIfile = ROIpath + "ROImap.zip";
				if( File.exists(ROIfile)){	
			           	print("Processing z-profiles for: " + input + file);
			           	//path = input + file;
			           	open(path);
			           	title=getTitle();
			           	      		
			           	
			           	
			           	roiManager("reset");
						
			           	//ROIpath = File.getParent(path) + File.separator;
						//ROIfile = ROIpath + "ROImap.zip";
						open(ROIfile);
						print("opened the ROIfile: " + ROIfile);
						ROIcount = roiManager("count");
						print("ROI count is :" + ROIcount);
						outerdir = File.getDirectory(input);
						tmpDir = "ROIoutputs" + File.separator;
						savefolder = File.getNameWithoutExtension(file)+ "_Results" + File.separator;
						if (File.exists(input + tmpDir) != 1) {
							//print("Making the ROIoutput directory.");
							File.makeDirectory(outerdir + tmpDir);}
							
						File.makeDirectory(outerdir + tmpDir + savefolder);
						tmpFile = File.getNameWithoutExtension(file);
			
					for(j=0;j<ROIcount;j++){
							selectWindow(title);
							roiManager("Select", j);
							run( "Plot Z-axis Profile" );
							Plot.getValues( x,y );
							lt = x.length;
							str = "";
							//print("Plotted Z-axis profile for ROI:" + j);
							
			           	
					for (i = 0; i < x.length; i++ ) {
						str += "" + x[i] + "," + y[i] + "\n"; }
					//print("Saving ROI.csv #" +j);	
					File.saveString( str, outerdir + tmpDir + savefolder + tmpFile   + "_" + "ROI" + j + ".csv" );	
					close();
					
					}
				cleanup();
				} else { 
					print("Didn't find an ROImap.zip");
					cleanup();
				}
				
	
    
    //} else {
    	//print("Found the wrong .tif file?");
    	//cleanup();
    //}
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


//extraneous functions from other steps in the data analysis pipeline below. 



// function to scan folders/subfolders/files to find files with correct suffix
function processFolder_ROIs(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder_ROIs("" + input + list[i]);
		if(endsWith(list[i], suffix))
			getROIs(input, list[i], dx);
	}
}

function getROIs(input, file, dx) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + file);
	path = input + file;
	open(path);
	title=getTitle();
	run("Find Maxima...", "prominence=3500 strict exclude output=[Single Points]");
			run("Analyze Particles...", "add"); 
			n = roiManager("count");
			if (n==0)
				exit("The ROI Manager is empty");
			for (i=0; i<n; i++) {
				roiManager("select", i);
				run("Enlarge...", "enlarge="+ dx + " pixel");
				roiManager("update");
			}
           	
           	
           	tmpDir = File.getDirectory(path) + File.separator;			
           	savefile = File.getNameWithoutExtension(path) + "ROImap.zip";
	
           	
           	output= tmpDir + savefile;
			roiManager("Save", output);
    		print("Saving ROIs to:" + output);   
	close(title);
	cleanup();
}




