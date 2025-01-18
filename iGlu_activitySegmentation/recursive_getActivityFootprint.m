%% This script combines elements from the CaImAn-MatLab package for processing and analyzing Ca2+ imaging videos 
%% with Philipe Mendonca's script for frequency-domain filtering of GluSnFR videos 
%% to achieve an automated script that loops through every .tiff in a folder and spits out max intensity projections.
%% A few caveats: 
%% memmap_file() was adopted as a solution to the memory problem in MatLab where the computer chokes on the massive videos I've produced from imaging spontaneous GluSnFR activity.
%% memmap_file() produces a data structure that represents the video file without requiring the video to sit in memory, thereby circumventing a key challenge. 
%% Still, memmap_file() doesn't play nice with .tif files > 4GB (a common issue) - after much Googling, it is currently above my paygrade to develop a script for that. Best practice is to avoid .tif >4GB altogether.


clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
clear;
format longg;
format compact;
%%% below, input your specific path for the MatLab scripts.
%addpath('C:\Users\blanpiedlab\Documents\MATLAB\ImarisReader');
addpath('\\blanpiedserver\NASShare3\Sam\Sam Matlab Scripts\activeScripts\iGlu_activitySegmentation');



%% setup path to file and package
%%% this script takes advantage of memmap_file() from CaImAn. This makes it
%%% easier to run the script, since you'll need to read your iGluSnFR
%%% videos into memory. The memmap file is a .mat file which will be
%%% created in the same folder as the .tif. The .mat file will be found on
%%% every subsequent run of the script if it exists, otherwise it will be
%%% generated de novo. 
%gcp;                                           % start local cluster
path_to_package = '\\blanpiedserver\NASShare3\Sam\Sam Matlab Scripts\activeScripts\iGlu_activitySegmentation';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
% Define the directory in which to find iGluSnFR videos. The script is recursive, so it will look in every subdirectory for .tif files.
start_path = fullfile('\\server\path\my_tif_files');
tic
% Ask user to confirm or change the directory for iGluSnFR videos.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end


% Get list of all subfolders.
allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {};
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);

for k = 1 : numberOfFolders
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	fprintf('Processing folder %s\n', thisFolder);
    % Get the TIFF files.
	filePattern = sprintf('%s/*.tif', thisFolder);
	baseFileNames = dir(filePattern);
    %print(filePattern);
	numberOfImageFiles = length(baseFileNames);
    disp(numberOfImageFiles);
    fprintf('testing the print function %s\n', filePattern);
    
    %% this part of the script is important. What is the framerate of your videos? frame_cycle is the time (in ms) each frame takes.
    frame_cycle=6.5; %given in ms; start with plate1, which was sampled at ET = 20 ms.
    period=frame_cycle/1000; %now in seconds    
    samplerate=1/period;
        
    mov_average=5; %Number of moving average (given in index, not time!)
	
if numberOfImageFiles >= 1
		% Go through all those image files.
        fprintf('Processing image files in folder %s\n', thisFolder)
		for f = 1 : numberOfImageFiles
			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
			fprintf('     Processing image file %s\n', fullFileName);
            
            imageStack_temp = [];
            is_memmaped = true;        % choose whether you want to load the file in memory or not

            %% load file

            if is_memmaped
                if exist([fullFileName(1:end-3),'mat'],'file')
                    data = matfile([fullFileName(1:end-3),'mat'],'Writable',true);
                else
                %sframe=150;						% user input: first frame to read (optional, default 1)
                %num2read=[];					% user input: how many frames to read   (optional, default until the end)
                %chunksize=[];                 % user input: read and map input in chunks (optional, default read all at once)
                data = memmap_file(fullFileName);
                %data = memmap_file_sequence(foldername);
                end
                sizY = size(data,'Y');                    % size of data matrix
            else
            T = 2000;                                 % load only a part of the file due to memory reasons
            data = read_file(fullFileName,1,T);
            sizY = size(data);
            end



            Y = double(data.Y);
            clearvars data
            im_filt=zeros(sizY(1),sizY(2),sizY(3));
                for i=1:size(Y,1)
                 %[f r i]% If you need to keep track of the processing time, uncomment this
                    for j=1:size(Y,2)
                        for k=1:size(Y,3)
                            y_tmp(k)=Y(i,j,k);
                        end
                    %Moving average: now that you have the trace (y_tmp) representing a single pixel in the z trajectory, filter it using a moving average.
                    y_smooth=movmean(y_tmp,mov_average);        
                    


                    % for JF646, 20 ms ET files, band-pass 0.35 - 200 Hz
                    % seems to work
                    % for GluSnFR3, 0.5 - 400 Hz? 
                    %Now that you have the trace (y_smooth) representing a single pixel in the z trajectory, filter it.
                    %% this is the important part of the script. Here, you should play with the numbers for the band-pass filter. 
                    %% For iGluSnFR videos with a 5 ms exposure time, I've found that low-pass filtering at 200 Hz and high-pass filtering at 0.5 Hz is effective. 
                    
                    
                    Ry_filt =y_smooth-mean(y_smooth(1:5)); %Subtract the baseline value (starting close to zero is the best way to avoid the filtering artifact)
                    Ry_filt =gaussfilter(Ry_filt,samplerate,200); % number is the low-pass filter i.e. <200 Hz
                    Ry_filt =Ry_filt-gaussfilter(Ry_filt,samplerate,0.1); % number is the high-pass filter i.e. > 0.1 Hz
                    
                    %add normalization step for binary maskings
                    %generate mean of each pixel
                    %generate standard deviation of each pixel
                    %(Ry-filt - mean_filt)/sd_filt
                    mean_filt = mean(y_smooth);
                    sd_filt = std(y_smooth);
                    full_filt = (Ry_filt - mean_filt)/sd_filt;
                    im_filt(i,j,:)=full_filt  ;
                    end 
                end
            clearvars imageStack_temp %just to optimse memory usage
            clearvars y_tmp y_smooth Y Ry_filt sizY sd_filt mean_filt full_filt
         
            %Scaling the corrected image
            im_filt(im_filt<0)=0;
            im_filt16=im_filt-min(im_filt(:));%Subtracting the minimum value, which will now be zero
            im_filt16=im_filt16.*(65535/max(im_filt16(:)));%Multiplying all values to 65535/max
     
            %For creating the corrected image
            im_filt16=uint16(im_filt16);
            im_filt16(im_filt16<1000)=0;%Offsets at 1000
 
            %Median filter
            med_filt='Yes'; condition=strcmp(med_filt,'Yes');
            if condition==1
                for i=1:size(im_filt16,3)
                    im_med16(:,:,i) = medfilt2(im_filt16(:,:,i));
                end
            else
            im_med16=im_filt16;    
         end
         
               
         %Now save the max projection of this stack
         im_med16_max=max(im_med16,[],3); %Now rescale to 16bit range       
         im_med16_max=uint16(im_med16_max);
         
         [folder, baseFileNameNoExt, extension] = fileparts(fullFileName);
         %saveImage = fullfile(topLevelFolder, saveDir);
         savePathFormat = [folder '\' baseFileNameNoExt '_maxFiltered.tiff'];
         
         imwrite(im_med16_max, savePathFormat);
         
         fprintf('     Finished processing image file %s\n', fullFileName);
         
         clearvars im_filt im_filt16 im_med16 im_med16_max
         end
	else
		fprintf('     Folder %s has no image files in it.\n', thisFolder);
	end
end
toc
