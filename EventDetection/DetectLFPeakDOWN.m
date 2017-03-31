function [ UPDOWNstates ] = DetectLFPeakDOWN( spiketimes,filtLFP,sf_LFP,detectionparms,detectints )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%detectionparms : minOFF,mininterOFF, peakthresh (STD), peakdist
%
%Current issue: does not account for spurious spikes in the middle of DOWN
%states...
%% DEV

%datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/WMData/';
%figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/AnalysisScripts/QW_DwellTimeMatchingAnalysis/';
%recname = 'SOMCell3';
%load(fullfile(datasetfolder,[recname,'.mat']))

%sf_LFP = 1./celldata(2).si;
%LFP = celldata(2).LFP;
%filtLFP = FiltNPhase(LFP,[1 4],sf_LFP);

%spiketimes = {celldata(2).MUA,celldata(2).spikes'};
%detectints = [celldata(2).NWh_on' celldata(2).NWh_off'];

t_LFP = (1:length(filtLFP))./sf_LFP;
minOFF = 0.025; %ms
mininterOFF = 0.03; %ms
peakthresh = 2.5; %for modZ
peakdist = 0.05;
%% Find peaks in the LFP


normLFP = NormToInt(filtLFP,detectints,sf_LFP,'modZ');

[peakheights,LFPeaks] = findpeaks(normLFP,t_LFP,'MinPeakHeight',peakthresh,'MinPeakDistance',peakdist);
[LFPeaks,keepPeaks] = RestrictInts(LFPeaks',detectints);LFPeaks=LFPeaks(:,1);
LFPeakheight = peakheights(keepPeaks);


%% Identify putative OFF states

allspikes = sort(cat(1,spiketimes{:}));

OFFints = [allspikes(1:end-1) allspikes(2:end)];
OFFdur = diff(OFFints,[],2);


%Find the OFFs around the peaks
[~,keeppeaks,keepOFF] = RestrictInts(LFPeaks,OFFints);
DOWNints = OFFints(keepOFF,:);
DOWNpeaks = LFPeaks(keeppeaks);
DOWNpeakmag = LFPeakheight(keeppeaks);
DOWNdur = diff(DOWNints,[],2);

%Remove DOWNs that are shorter than the threshold
removepeaks = DOWNdur<minOFF;
DOWNints(removepeaks,:) = [];DOWNdur(removepeaks) = [];
DOWNpeaks(removepeaks) = []; DOWNpeakmag(removepeaks) = [];

%Calculate UPs and remove those that aren't in detectionints
UPints = [DOWNints(1:end-1,2) DOWNints(2:end,1)];
UPints = RestrictInts(UPints,detectints);

UPdur = diff(UPints,[],2);

%% Ouput in .event.mat format

detectionparms.minOFF = minOFF; %ms
detectionparms.mininterOFF = mininterOFF; %ms
detectionparms.peakthresh = peakthresh; %for modZ
detectionparms.peakdist = peakdist;


UPDOWNstates.ints.UP = UPints;
UPDOWNstates.ints.DOWN = DOWNints;
UPDOWNstates.DOWNpeaks = DOWNpeaks;
UPDOWNstates.DOWNpeakmag = DOWNpeakmag;
UPDOWNstates.detectorinfo.detectorname = 'DetectLFPeakDOWN';
UPDOWNstates.detectorinfo.detectionparms = detectionparms;


%%
% figure
% subplot(2,2,1)
%     hist(log10(UPdur),50)
%     LogScale('x',10)
% subplot(2,2,2)
%     hist(log10(DOWNdur),50)
%     LogScale('x',10)


end

