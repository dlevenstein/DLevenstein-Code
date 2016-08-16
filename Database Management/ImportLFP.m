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
display(['Importing: ',num2str(ss), ' of ', num2str(numrecs)])
datasetrecfolder = fullfile(datasetfolder{ss},recordingname{ss});
sourcerecfolder = fullfile(sourcefolder{ss},recordingname{ss});
if ~exist(datasetrecfolder,'dir')
    mkdir(datasetrecfolder)
end

    %GoodSleepInterval for BWData
    goodsleepmatBW = fullfile(sourcerecfolder,[recordingname{ss},'_GoodSleepInterval.mat']);
    goodsleepmatDL = fullfile(datasetrecfolder,[recordingname{ss},'_GoodInterval.mat']);
    if exist(goodsleepmatBW,'file')
        load(goodsleepmatBW)
      	scoretime = StartEnd(GoodSleepInterval,'s');
        copyfile(goodsleepmatBW,datasetrecfolder)
    elseif exist(goodsleepmatDL,'file')
        load(goodsleepmatDL)
      	scoretime = GoodSleepInterval;
    else
        scoretime = [0 Inf];
    end
    
    %Bad Channels - copy to source folder.  Need to deal with this for
    %version control.......
    badchannelDL = fullfile(datasetrecfolder,'bad_channels.txt');
    if exist(badchannelDL,'file') 
        copyfile(badchannelDL,sourcerecfolder,'f')
    end
    
    %SleepScore the data from source and save in dropbox database
    SleepScoreMaster(sourcefolder{ss},recordingname{ss},...
        'savedir',datasetfolder{ss},'spindledelta',false,...
        'scoretime',scoretime,'overwrite',true)
    close all
    
    
    %Add regional LFP to the database using SpikeGroupAnatomy, remove 60Hz noise
    SGAfilename_source = fullfile(sourcerecfolder,[recordingname{ss},'_SpikeGroupAnatomy.csv']);
    SGCfilename_dataset = fullfile(datasetrecfolder,[recordingname{ss},'_SpikeGroupAnatomy.csv']);
    if exist(SGAfilename_source,'file') & ~exist(SGCfilename_dataset,'file')
        copyfile(SGAfilename_source,datasetrecfolder)
    elseif ~exist(SGAfilename_source,'file') & ~exist(SGCfilename_dataset,'file')
        display('No Spikegroup Anatomy File...')
        continue
    end
    
    scoremetricspath = fullfile(datasetrecfolder,[recordingname{ss},'_StateScoreMetrics.mat']);
    load(scoremetricspath,'SWchannum','THchannum')
    
    CleanAndCopyRegionalLFP(datasetfolder{ss},recordingname{ss},sourcefolder{ss},...
        [SWchannum,THchannum],scoretime)

end


end

