function [  ] = CleanAndCopyRegionalLFP(datasetfolder,recname, sourcefolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%

xmlfilename = fullfile(sourcefolder,recname,[recname,'.xml']);
if exist (fullfile(sourcefolder,recname,[recname,'.lfp']),'file')
    rawlfppath = fullfile(sourcefolder,recname,[recname,'.lfp']);
elseif exist (fullfile(sourcefolder,recname,[recname,'.eeg']),'file')
    rawlfppath = fullfile(sourcefolder,recname,[recname,'.eeg']);
else 
    display('No .lfp file')
end


Par = LoadPar_SleepScore(xmlfilename);
Fs = Par.lfpSampleRate; % Hz, LFP sampling rate
nChannels = Par.nChannels;

if isfield(Par,'SpkGrps')
    SpkGrps = Par.SpkGrps;
elseif isfield(Par,'AnatGrps')
    SpkGrps = Par.AnatGrps;
    display('No SpikeGroups, Using AnatomyGroups')
else
    display('No SpikeGroups...')
end

%% Load the Region Name List
%Load region name list
laptopdatabasefolder = '/Users/dlevenstein/Dropbox/Research/Datasets';
desktopdatabasefolder = '/mnt/data1/Dropbox/research/Datasets';
if exist(desktopdatabasefolder,'dir')
    dropboxdatabasefolder = desktopdatabasefolder;
elseif exist(laptopdatabasefolder,'dir')
     dropboxdatabasefolder = laptopdatabasefolder;
else
    display('No Dropbox Folder...')
end
regionnamefilename = fullfile(dropboxdatabasefolder,'RegionNameList.csv');
regionnames=readtable(regionnamefilename);
numregions = width(regionnames);

%% Load the channel anatomy map
spikegroupanatomyfilename = fullfile(datasetfolder,recname,[recname,'_SpikeGroupAnatomy.csv']);
if ~exist(spikegroupanatomyfilename,'file')
    display('No SpikeGroup Anatomy file...')
    return
end
spkgroupanatomy=readtable(spikegroupanatomyfilename);
%Find unique site identifiers.
sitenames = unique(spkgroupanatomy.AnatomicalSite);

for rr = 1:numregions
    regionlabels = regionnames.(rr);
    regionroups = ismember(sitenames,regionlabels);
    %-find the sitenames that fit within the region
    %-for each of those sites, find the channel numbers that correspond to
    %that spikegroup, exclude bad_channels
    %
end

ctxgroups = ismember(spkgroupanatomy.AnatomicalSite,CTXlabels);
ctxgroups = spkgroupanatomy.SpikeGroup(ctxgroups);
ctxchannels = [spikegroups{ctxgroups}];

end

