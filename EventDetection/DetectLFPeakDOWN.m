function [ UPDOWNstates ] = DetectLFPeakDOWN( basePath,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   basePath    
%   'NREMInts'  Interval of times for NREM 
%               (Default: loaded from SleepState.states.mat)
%   'SWChann'   Channel with the most robust (positively deflecting) Slow
%               Waves. (0-Indexing a la neuroscope)
%               (Default: use SWChannel from sleep scoring)
%   'CTXSpkGroups' Spike groups that are in the cortex...  default: all
%   'detectionparms' Structure of detection parameters
%
%detectionparms : minOFF,mininterOFF, peakthresh (STD), peakdist (currently
%hardcoded)
%
%Current issue: does not account for spurious spikes in the middle of DOWN
%states...

%% DEV for buzcode...
% baseName = 'DT2_rPPC_rCCG_362um_218um_20160209_160209_183610';
% baseFolder ='/mnt/proraidDL/Database/DTData/';
% LFPchan = 12;
% SpikeGroups = 1:5
% 
% 
% lfpName = fullfile(baseFolder,baseName,[baseName.lfp]);
% 

%Defaults
basePath = pwd;
CTXSpikeGroups = 'all';
SWChann = [];
NREMInts = [];

%Parms
minOFF = 0.025; %ms
mininterOFF = 0.03; %ms
peakthresh = 2.5; %for modZ
lowerpeakthresh = 1; %for modZ
peakdist = 0.05;


%% Collect all the Necessary Pieces of Information

%Sleep Scoring (for NREM)
[SleepState] = bz_LoadStates(basePath,'SleepState');
if isempty(SleepState) && isempty(NREMInts)
    button = questdlg(['SleepState.states.mat does not exist, '...
        'would you like to run SleepScoreMaster?'],...
        'DetectLFPeakDOWN Needs NREM',...
        'Yes','No, use all time','Cancel','Yes');
    switch button
        case 'Yes'
            SleepState = SleepScoreMaster(basePath);
            NREMInts = SleepState.ints.NREMstate;
        case 'No, use all time'
            NREMInts = [0 Inf];
            SleepState.detectorparams.empty = [];
        case 'Cancel'
            return
    end
else
    NREMInts = SleepState.ints.NREMstate;
end
   
%LFP in the detection channel
if ~isfield(SleepState.detectorparams,'SWchannum') && isempty(SWChann)
    SWChann = input('Which channel shows the most robust slow waves?');
elseif isempty(SWChann)
    SWChann = SleepState.detectorparams.SWchannum;
end
lfp = bz_GetLFP(SWChann,'basepath',basePath);

%Spikes in the CTX Spike Groups
spikes = bz_GetSpikes('basepath',basePath,'spikeGroups',CTXSpikeGroups);

%% Filter the LFP

filterbounds = [0.5 4];
deltaLFP = bz_Filter(lfp,'passband',filterbounds,'filter','fir1','order',3);

%%
%deltaLFP = FiltNPhase(LFP,deltarange,sf,numcyc);


%% Find peaks in the LFP


normLFP = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

[peakheights,LFPeaks] = findpeaks(normLFP,deltaLFP.timestamps,'MinPeakHeight',lowerpeakthresh,'MinPeakDistance',peakdist);
[LFPeaks,keepPeaks] = RestrictInts(LFPeaks,NREMInts);LFPeaks=LFPeaks(:,1);
LFPeakheight = peakheights(keepPeaks);

%%
numcells = length(spikes.times);
allspikes = sort(cat(1,spikes.times{:}));
allspikes = sort(cat(1,S_CellFormat{:})); %For .mat... update once fixed.
numcells = length(S_CellFormat);
%%
win = 1;
reltime = [];
spkpeakheight = [];
for pp = 1:length(LFPeaks)
    pp
    nearpeakspikes = allspikes >= LFPeaks(pp)-win & allspikes <= LFPeaks(pp)+win;
    reltime = [reltime; allspikes(nearpeakspikes)-LFPeaks(pp)];
    spkpeakheight = [spkpeakheight; LFPeakheight(pp).*ones(sum(nearpeakspikes),1)];
   % keyboard
end

%% Spike-SW CCG histogram by peak height
numtimebins = 200;
nummagbins = 30;
timebins = linspace(-win,win,numtimebins);
magbins = linspace(1,8, nummagbins);
peakmagdist = hist(LFPeakheight,magbins);

spikehitmat = hist3([reltime,spkpeakheight],{timebins,magbins});
ratemat = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);

%% Peak threshold for putative SWs

putSWs = LFPeaks(LFPeakheight>=peakthresh);
putSWPeakHeights = LFPeakheight(LFPeakheight>=peakthresh);

%% Calculate Binned Population Rate
dt = 0.003; %dt = 3ms
overlap = 10; %Overlap = 10 dt
winsize = dt*overlap; %meaning windows are 30ms big
[spikemat,t,spindices] = SpktToSpkmat(allspikes, [], dt,overlap);
 SpktToSpkmat

%% CCG of SLow Waves
[SW_CCG,t_CCG] = CCG(putSWs,ones(size(putSWs)),'binSize',0.02);

%%
figure
bar(t_CCG,SW_CCG)
%%
ratecolor = makeColorMap([1 1 1],[0.8 0 0],[0 0 0]);
figure
subplot(2,2,1)
    plot(reltime,(spkpeakheight),'k.','markersize',1)
subplot(2,2,4)
    imagesc(timebins,magbins,ratemat')
    hold on
    plot(get(gca,'xlim'),peakthresh.*[1 1],'k--')
    plot([0 0],get(gca,'ylim'),'k-')
    colormap(ratecolor)
    colorbar('location','east')
    xlim([-0.8 0.8])
    axis xy
    xlabel('t (relative to delta peak)');ylabel('Delta Peak Amplitude (modZ)')
    
subplot(4,2,2)
    bar(t_CCG,SW_CCG)
    xlim([-0.8 0.8])
%%
winsize = 10; %s
samplewin = randsample(LFPeaks,1)+winsize.*[-0.5 0.5];

figure
subplot(4,1,1)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    box off
    plot(lfp.timestamps(ismember(lfp.timestamps,LFPeaks)),lfp.data(ismember(lfp.timestamps,LFPeaks)),'r.')
    xlim(samplewin)
subplot(4,1,2)
    plot(deltaLFP.timestamps,normLFP,'b')
    hold on
    xlim(samplewin)
    plot(get(gca,'xlim'),peakthresh.*[1 1],'r--')
%plot(deltaLFP.timestamps,deltaLFP.data,'g')
%plot(deltaLFP.timestamps(ismember(deltaLFP.timestamps,LFPeaks)),normLFP(ismember(deltaLFP.timestamps,LFPeaks)),'r.')


%% Identify putative OFF states

LFPeaks = LFPeaks(LFPeakheight>=peakthresh);
LFPeakheight = LFPeakheight(LFPeakheight>=peakthresh);
%allspikes = sort(cat(1,spiketimes{:}));

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
UPints = RestrictInts(UPints,NREMInts);

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

