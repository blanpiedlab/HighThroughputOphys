//Script structure will process 1 folder full of ims files using the subfunction "processfiles"
//Macro calls in processfiles can be edited to whatever you need done, just use the macro recorder
//
//As written, this will convert a folder of .ims files (timelapses) to .tiff files 
//NOTE: This script won't work together with MatLab unless the image files are sufficiently small! <4 GB

//Select input folders
input = getDirectory("Please choose the folder with your registered and cropped .tif files"); 
input = replace(input,"/",File.separator); //replace any / with \ in path name
prefix="area";
suffix=".tif"
//This returns the path of the parent folder with slashes in this direction \

//Batchmode makes it so that the images don't open, which is usually faster. 
//If you want to see the images open and process then change setBatchMode to false.
//Sometimes something you have coded with macro recorder won't work in batch mode...


//The for loop will iterate through all items in the input directory.
setBatchMode(true);
processFolder_splitFiles(input)
print("Analysis is done");
beep();
setBatchMode(false);





function processFolder_splitFiles(input) {
	list = getFileList(input);
	for (i=0; i<list.length; i++) {
    	if(File.isDirectory(input + list[i]))
			processFolder_splitFiles("" + input + list[i]);
    	if (startsWith(list[i], prefix) && endsWith(list[i], suffix)) {//loop will find all tifs
    		full_filename = list[i];
    		
    		processfiles(input,full_filename); //runs function processfiles on the file //deleted c3 from function call
        
		}
	}
}




//Edit this function as needed for your workflow
function processfiles(input,filename) { 
	print("Processing file: " + filename);
	//Import the highest resolution (series_1) file as a hyperstack with composite channel settings
	open(input + filename);
	origIm = getTitle();
	selectWindow(origIm);
	run("Duplicate...", "duplicate channels=1");
	GluIm = "C1-"+File.getNameWithoutExtension(origIm) + ".tif";
	saveAs("Tiff", input + GluIm); //save the tif in the new folder. you can edit this to change filenames if desired
	
	selectWindow(origIm);
	run("Duplicate...", "duplicate channels=2");
	JFIm = "C2-"+File.getNameWithoutExtension(origIm) + ".tif";
	saveAs("Tiff", input +  JFIm); //save the tif in the new folder. you can edit this to change filenames if desired
	
	
	
	//Save file
	
	run("Close All"); //close all windows
	
}