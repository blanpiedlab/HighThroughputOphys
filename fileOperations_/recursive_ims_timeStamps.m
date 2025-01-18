% Start with a folder and get a list of all subfolders.
% Finds and prints names of all PNG, JPG, and TIF images in 
% that folder and all of its subfolders.
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
format longg;
format compact;
addpath('\\blanpiedserver\NASShare3\Sam\Sam Matlab Scripts\github_repo');
% 
% Define a starting folder.
start_path = fullfile('\\blanpiedserver\NASShare3\Sam\2021-11-15_toBeAnalyzed');
tic
% Ask user to confirm or change.
topLevelFolder = uigetdir(start_path);
if topLevelFolder == 0
	return;
end

%make a directory in which to store .csv output
saveDir = 'timeStamps';
mkdir(topLevelFolder, saveDir);

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
numberOfFolders = length(listOfFolderNames)

for k = 1 : numberOfFolders
	% Get this folder and print it out.
	thisFolder = listOfFolderNames{k};
	fprintf('Processing folder %s\n', thisFolder);
    % Get the IMS files.
	filePattern = sprintf('%s/*.ims', thisFolder);
	baseFileNames = dir(filePattern);
	numberOfImageFiles = length(baseFileNames);
	% Now we have a list of all files in this folder.
	
	if numberOfImageFiles >= 1
		% Go through all those image files.
		for f = 1 : numberOfImageFiles
			fullFileName = fullfile(thisFolder, baseFileNames(f).name);
			fprintf('     Processing image file %s\n', fullFileName);
            
            %do ImarisReader things
            imObj = ImarisReader(fullFileName);
	
            timeStamps = imObj.DataSet.Timestamps;
    
            %convert timeStamps to time 0 start 
            t = datenum(timeStamps,'yyyy-mm-dd HH:MM:SS.FFF');
            dt = (t - t(1))*24*3600; %convert to seconds
            colNames = {'timeStamps','dt'};
    
            imTable = table(timeStamps, dt, 'VariableNames', colNames);
    
            %get the filename
            [~,fid,~] = fileparts(fullFileName);    
            %get identifiers from the filename and add to the table
            imTable.fileID(1:height(imTable)) = {fid};
            imTable.index = (1:height(imTable)).' ;
    
            %reorder the table
            imTable = imTable(:,{ 'fileID' 'index' 'timeStamps' 'dt'});
            %saving the table
            [folder, baseFileNameNoExt, extension] = fileparts(fullFileName);
            saveTable = fullfile(topLevelFolder, saveDir);
            savePathFormat = [saveTable '\' baseFileNameNoExt '_timestamps.csv'];
            writetable(imTable,savePathFormat,'Delimiter',',');
            
            
		end
	else
		fprintf('     Folder %s has no image files in it.\n', thisFolder);
	end
end
toc