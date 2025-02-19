//Script structure will process 1 folder full of ims files using the subfunction "processfiles"
//Macro calls in processfiles can be edited to whatever you need done, just use the macro recorder
//
//As written, this will convert a folder of .ims files (timelapses) to .tiff files 
//NOTE: This script won't work together with MatLab unless the image files are sufficiently small! <4 GB

//Select input folders
input = getDirectory("Please choose the folder with your ims files"); 
input = replace(input,"/",File.separator); //replace any / with \ in path name

//This returns the path of the parent folder with slashes in this direction \

//Batchmode makes it so that the images don't open, which is usually faster. 
//If you want to see the images open and process then change setBatchMode to false.
//Sometimes something you have coded with macro recorder won't work in batch mode...


//The for loop will iterate through all items in the input directory.
setBatchMode(true);
list = getFileList(input); 
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".ims")) { //loop will find files that end with .ims
    	count=i;
    	filename = list[i];
    	print("Processing " + filename);
    	
        processfiles(input,filename,count); //runs function processfiles on the file
	}
}
print("Analysis is done");
beep();
setBatchMode(false);

//Edit this function as needed for your workflow
function processfiles(input,filename, count) { 
	
	//Import the highest resolution (series_1) file as a hyperstack with composite channel settings
	run("Bio-Formats Importer", "open=[" + input + filename + "] color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	//Rename the image to just the base filename rather than the long thing that ims series loads as
	//Otherwise you will have issues with filenames getting cut off
	x = indexOf(filename,".ims");
	origIm = substring(filename,0,x);
	origIm = replace(origIm," ","_");
	rename(origIm);
	n= nSlices;
	last = n; //n-1
	
	//run("Slice Keeper", "first=1 last="+last+" increment=1"); //added functionality to delete final frame so the weird mirror flip on last frame in my videos from 12.14.21 doesn't impact analysis pipeline
	run("Subtract Background...","rolling=50 stack"); //add rolling ball subtract background prior to GaussFilter
	
	savefolder = "_Results_" + count + File.separator;
	File.makeDirectory(input + savefolder);


	//Save file
	saveAs("Tiff", input + savefolder + origIm); //save the tif in the new folder. you can edit this to change filenames if desired
	run("Close All"); //close all windows
	
}
