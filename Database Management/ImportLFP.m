function [ ] = ImportLFP()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%

%%  Select Recordings from DatasetGuide spreadsheet
laptopdatabasefolder = '/Users/dlevenstein/Dropbox/Research/Datasets';
desktopdatabasefolder = '/mnt/data1/Dropbox/research/Datasets';
if exist(desktopdatabasefolder,'dir')
    dropboxdatabasefolder = desktopdatabasefolder;
elseif exist(laptopdatabasefolder,'dir')
     dropboxdatabasefolder = laptopdatabasefolder;
else
    display('No Dropbox Folder...')
end

datasetguidefilename = fullfile(dropboxdatabasefolder,'DatasetGuide.xlsx');
datasetguide=readtable(datasetguidefilename);
possiblerecordingnames = datasetguide.RecordingName;
possibledatasetfolders = datasetguide.Dataset;
[s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                'ListString',possiblerecordingnames);
recordingname = possiblerecordingnames(s);
datasetfolder = cellfun(@(X) fullfile(dropboxdatabasefolder,X),...
    possibledatasetfolders(s),'UniformOutput',false);
sourcefolder = datasetguide.DesktopFolder(s);


numrecs = length(s);
%% Sleep Score the LFP from the source folder

for ss = 1:numrecs
    
    %GoodSleepInterval for BWData
    goodsleepmatBW = fullfile(sourcefolder{ss},recordingname{ss},[recordingname{ss},'_GoodSleepInterval.mat']);
    if exist(goodsleepmatBW,'file')
        load(goodsleepmatBW)
        
      %  scoretime = 
    else
        scoretime = [0 Inf];
    end
    
    %SleepScore the data from source and save in dropbox database
    SleepScoreMaster(sourcefolder{ss},recordingname{ss},...
        'savedir',datasetfolder{ss},'spindledelta',true,...
        'scoretime',scoretime)
    close all
end


end

