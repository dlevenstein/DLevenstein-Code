function [  ] = CleanAndCopyRegionalLFP(datasetfolder,recname,sourcefolder,specialchannels,scoretime)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
PLOTFIG = true;
figstuff = struct;
figfolder = fullfile(datasetfolder,recname);

if ~exist('sourcefolder','var')
    sourcefolder = datasetfolder;
end

xmlfilename = fullfile(sourcefolder,recname,[recname,'.xml']);
if exist (fullfile(sourcefolder,recname,[recname,'.lfp']),'file')
    rawlfppath = fullfile(sourcefolder,recname,[recname,'.lfp']);
elseif exist (fullfile(sourcefolder,recname,[recname,'.eeg']),'file')
    rawlfppath = fullfile(sourcefolder,recname,[recname,'.eeg']);
else 
    display('No .lfp file')
end

lfpsavename = fullfile(datasetfolder,recname,[recname,'_LFP.mat']);

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

%Load Rejectchannels list
rejectchannels = [];
if exist(fullfile(datasetfolder,recname,'bad_channels.txt'),'file')%bad channels is an ascii/text file where all lines below the last blank line are assumed to each have a single entry of a number of a bad channel (base 0)
    t = ReadBadChannels_SleepScore(fullfile(datasetfolder,recname));
    rejectchannels = cat(1,rejectchannels(:),t(:));
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


%% Load the channel anatomy map
spikegroupanatomyfilename = fullfile(datasetfolder,recname,[recname,'_SpikeGroupAnatomy.csv']);
if ~exist(spikegroupanatomyfilename,'file')
    display('No SpikeGroup Anatomy file...')
    return
end
spkgroupanatomy=readtable(spikegroupanatomyfilename);
%Find unique site identifiers.
[sitenames,~,siteidx] = unique(spkgroupanatomy.AnatomicalSite);
%Find and site names that aren't in a region list... need to save these as
%separate region
allsubregionnames = {regionnames{:,:}{:}};
allsubregionnames(strcmp('',allsubregionnames))=[]; %remove empty strings
notonlist = sitenames(~ismember(sitenames,allsubregionnames));
notonlist(strcmp('-',notonlist))=[]; %remove don't-count channels
notonlist(strcmp('',notonlist))=[]; %remove empty strings
for nn = 1:length(notonlist)
    regionnames.(notonlist{nn}) = cell(size(regionnames.(1)));
    regionnames.(notonlist{nn})(1) = notonlist(nn);
    regionnames.(notonlist{nn})(2:end)={''};
end
%%

%%
numregions = width(regionnames);
for rr = 1:numregions
%%
%rr = 1;
    regionlabels = regionnames.(rr);
    regiongroups = find(ismember(sitenames,regionlabels));
    numinregiongroups = length(regiongroups);
    %%
    for gg = 1:numinregiongroups
      %  gg = 1;
        groupname = sitenames{regiongroups(gg)};
        ingroupspikegroups = find(siteidx==regiongroups(gg));
        groupchannels = [SpkGrps(ingroupspikegroups).Channels];
        
        %for each of those sites, find the channel numbers that correspond to
        %that spikegroup, exclude bad_channels, select usechannels
        if sum(ismember(specialchannels,groupchannels))>0
            groupchannels = intersect(specialchannels,groupchannels);
            if ismember(groupchannels,rejectchannels)
                display('special channel is marked as a bad channel, you might want to check this...')
            end
        else
            groupchannels = setdiff(groupchannels,rejectchannels);
        end
    
    %%
        %Load the right LFP
        pLFP = LoadBinary_Down(rawlfppath,'frequency',Fs,...
        'nchannels',nChannels,'channels',groupchannels+1,...
        'start',scoretime(1),'duration',diff(scoretime));
        %%
        %Compare LFP channels
        [noisescore,noiseresults] = CompareLFPChannels(pLFP);
        [~,minnoisechannel] = min(noisescore);
        pLFP = pLFP(:,minnoisechannel);
        noisetimes = noiseresults.noisetimes(:,minnoisechannel);
        %Pick the channel with the least noise (other criteria as well?
        %include in CompareLFPChannels
    %%
        %-Remove 60Hz Noise

        %add LFP to LFP structure
        regionname = regionnames.Properties.VariableNames{rr};
        if gg ~= 1
            regionname = [regionname,num2str(gg)];
        end
        LFP.(regionname) = pLFP;
        LFP.subregionname.(regionname) = groupname;
        LFP.channum.(regionname) = groupchannels(minnoisechannel);
        LFP.noisetimes.(regionname) = find(noisetimes);
        LFP.spikegroups.(regionname) = ingroupspikegroups;
        
        
        %Stuff for plotting
        if PLOTFIG
            downfactor = 5;
            figstuff(end+1).regionname = regionname;
            figstuff(end).subregionname = groupname;
            figstuff(end).LFP = zscore(downsample(pLFP,downfactor));
            %Calcualte Z-scored Spectrogram
            numfreqs = 100;
            freqlist = logspace(0,2,numfreqs);
            window = 10;
            noverlap = 9;
            window = window*Fs;
            noverlap = noverlap*Fs;
            [FFTspec,FFTfreqs,t_FFT] = spectrogram(figstuff(end).LFP,window,noverlap,freqlist,Fs./downfactor);
            figstuff(end).FFTspec = log10(abs(FFTspec));
            [~,figstuff(end).mu,figstuff(end).sig] = zscore((figstuff(end).FFTspec)');
        end
    end

end
figstuff(1) = [];

LFP.sf = Fs;
LFP.t = [1:length(pLFP)]./LFP.sf;
%%
save(lfpsavename,'LFP')


%% Figure
if PLOTFIG
viewwin  =[t_FFT(1) t_FFT(end)];
t_LFP = downsample(LFP.t,downfactor);
    fig = figure;
    suptitle([recname,': Selected Regional LFP'])
    for ff = 1:length(figstuff)
        subplot(3,1,ff)
            imagesc(t_FFT,log2(FFTfreqs),(figstuff(ff).FFTspec))
            hold on
            plot(t_LFP,((figstuff(ff).LFP)./10)-1,'k')
            axis xy
            set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
            set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
            caxis([3.5 6.5])
            caxis([min(figstuff(ff).mu)-2.5*max(figstuff(ff).sig),...
                max(figstuff(ff).mu)+2.5*max(figstuff(ff).sig)])
            xlim(viewwin)
           % colorbar('east')
            ylim([log2(FFTfreqs(1))-2 log2(FFTfreqs(end))])
            box off
            set(gca,'XTickLabel',{})
            ylabel({[figstuff(ff).regionname,': ',figstuff(ff).subregionname],'f (Hz)'})
            
            if mod(ff,3) == 0 | ff==length(figstuff)
                 saveas(fig,[figfolder,'/',recname,'_RegionLFP'],'jpeg')
                 figure
            end
    end
    close all
end

end

