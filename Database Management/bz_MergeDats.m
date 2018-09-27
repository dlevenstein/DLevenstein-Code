function [  ] = bz_CatDats( topfolder )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% For now, only works in linux
%
%
%Tries to merge the following .dat files: 
%amplifier, supply, time, auxiliary, digitalin, analogin
%
%If no amplifier.dat file is found, tries baseName.dat
%
%If no .xml is found in the top directory, copies from one of the
%subfolders
%
%
%%
topfolder = '/mnt/proraidDL/Database/WMProbeData/180211_WT_M1M3_LFP_Layers_Pupil_EMG_Pole';
baseName = bz_BasenameFromBasepath(topfolder);
%%
[basePaths, dirNames, fileNames] = dirwalk(topfolder);
fullfolders = cellfun(@(X) [basePaths{1},filesep,X],dirNames{1},'uniformoutput',false);
folders = dirNames{1};
%allfiles = fileNames(ismember(basePaths,folders));
datfiles = cellfun(@(X) dir([X,filesep,'*.dat']),fullfolders,'UniformOutput',false);
datfiles = cellfun(@(X) {X(:).name},datfiles,'UniformOutput',false);
%%
dattypes = {'amplifier','supply','time','auxiliary','digitalin','analogin'}';

for dd = length(dattypes):-1:1
    whichdats = cellfun(@(X) ismember(X,[dattypes{dd},'.dat']),datfiles,'UniformOutput',false);
    numtypefiles.(dattypes{dd}) = cellfun(@(X) sum(X),whichdats);
    
    if all(numtypefiles.(dattypes{dd})==0)
        display(['No ',dattypes{dd},' files. Skipping...'])
        dattypes(dd) = [];
    elseif any(numtypefiles.(dattypes{dd})==0)
        display(['Some of your folders have no ',dattypes{dd},' file. Skipping those folders...'])
        %Skip those folders
    else
        fullfilenames.(dattypes{dd}) =  cellfun(@(X,Y,Z) fullfile(X,Y(Z)),fullfolders,datfiles,whichdats,'UniformOutput',false);
        localfilenames.(dattypes{dd}) = cellfun(@(X,Y,Z) fullfile(X,Y(Z)),folders,datfiles,whichdats,'UniformOutput',false);
    end
end

%% Load the xml
try 
    %Look for xml/sessionInfo in topfolder
    sessionInfo = bz_getSessionInfo(topfolder,'noPrompts',true);
catch
    %If none exists, look for xml in any of the subpaths
    disp('No .xml or .sessionInfo in top folder, trying subfolders')
    for ff = 1:length(fullfolders)
        try
            [sessionInfo,xmlfilename] = LoadParameters(fullfolders{ff});
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(...
                xmlfilename,...
                fullfile(topfolder,[baseName,'.xml']),'f');
            display(['Copied xml from ',folders{ff}])
            break
        catch
        end
    end
end

%% Get time points from the time.dat
%Use the timestamps from time.dat to get the sort order
%Number of samples in time.dat. First timepoint, last timepoint
%
%Convert from number of samples to recording time of start/ends
% 
for ff = 1:length(fullfolders)
%ff = 1;
	timefilename = fullfile(fullfolders{ff},'time.dat');
	f = fopen(timefilename,'r'); 
    % Determine total number of samples in file
    fileStart = ftell(f);
    
    %Read the first time point
    firsttimepoint = fread(f,1,'int32');
    status = fseek(f,-4,'eof'); %int32 = 4 bytes
    lasttimepoint = fread(f,1,'int32');
%     if status ~= 0,
%         fclose(f);
%         error('Error reading the data file (possible reasons include trying to read past the end of the file).');
%     end
    fileStop = ftell(f);
    
    firstlasttimepoints(ff,:) = [firsttimepoint lasttimepoint];
    numsamples(ff) = fileStop./4;
    if ff==1
        transitiontimes = firstlasttimepoints;
    else
        transitiontimes(ff,:) = firstlasttimepoints(ff,:)+firstlasttimepoints(ff-1,2)+1;
    end
end

transitiontimes_s = transitiontimes./sessionInfo.rates.wideband; %convert to seconds


%% Merge the .dats!




%% Events file
eventsfilename = fullfile(topfolder,[baseName,'.MergePoints.events.mat']);

MergePoints.timestamps = transitiontimes_s;
MergePoints.timestamps_sampels = transitiontimes;
MergePoints.foldernames = folders;
MergePoints.filesmerged = localfilenames;
MergePoints.detectorinfo.detectorname = 'bz_MergeDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

%Saving SleepStates
save(eventsfilename,'MergePoints');
    %%
%end


%%
% supplyfiles = dir([topfolder,'/*/supply.dat']);
% timefiles = dir([topfolder,'/*/time.dat']);
% auxiliaryfiles = dir([topfolder,'/*/auxiliary.dat']);
% digitalinfiles = dir([topfolder,'/*/digitalin.dat']);
% analoginfiles = dir([topfolder,'/*/analogin.dat']);


%Merge the dats
%%

%unix('*/supply.dat')
%Make an events.mat file indicating the transition times (clock and
%recording time)
%Write the transition times and other things to an events.mat file
%% Check the time.dat file for switch times and that the
% dat files have the right number of samples to match the time file?
end

