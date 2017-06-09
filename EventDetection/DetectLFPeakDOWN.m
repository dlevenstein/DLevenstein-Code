function [ UPDOWNstates ] = DetectLFPeakDOWN( basePath,varargin)
%[UPDOWNstates] = DetectLFPeakDOWN(basePath) detects neocortical slow
%waves using a combination of a positive deflection in the LFP (delta wave)
%and period of neuronal inactivity (DOWN/OFF state).
%
%INPUTS
%   basePath    
%   'NREMInts'  -Interval of times for NREM 
%               -(Default: loaded from SleepState.states.mat)
%   'SWChann'   -Channel with the most robust (positively deflecting) Slow
%               Waves. (0-Indexing a la neuroscope)
%               (Default: use SWChannel from sleep scoring)
%   'CTXSpkGroups' -Spike groups that are in the cortex...  default: all
%   'detectionparms' -Structure of detection parameters
%       .minOFF         (default: 0.025s)
%       .mininterOFF    (default: 0.03s)
%       .peakdist       (default: 0.05s)
%       .peakthresh     (default: 2.5 stds - modified Z score)
%   'SHOWFIG'   -true/false show a quality control figure (default: true)
%   'saveMat'   -logical (default=true) to save in buzcode format
%
%detectionparms : minOFF,mininterOFF, peakthresh (STD), peakdist (currently
%hardcoded)
%
%OUTPUTS
%   UPDOWNstates    a buzcode structure
%
%
%
%DLevenstein 2016/2017
%TO DO
%-finish commenting and input parsing...
%-incorporate multiple channels for detection of slow wave, which is robust
%on all (deep) lfp channels in the local cortical population
%-incorporate the drop in gamma power (a la Watson et al 2016). especially
%useful for recordings with low cell count
%% Defaults and Parms


%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end
CTXSpikeGroups = 'all';
SWChann = [];
NREMInts = [];
SAVEMAT = true;
FORCERELOAD = false;
SHOWFIG = true;

%Parms
minOFF = 0.025; %ms
mininterOFF = 0.03; %ms
peakthresh = 2.5; %for modZ
lowerpeakthresh = 1; %for modZ
peakdist = 0.05;

goodcellnumber = 20;

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.UPDOWNstates.states.mat']);

if exist(savefile,'file') && ~FORCERELOAD
    display(['Slow Oscillation alter Detected, loading ',baseName,'.UPDOWNstates.states.mat'])
    UPDOWNstates = bz_LoadStates(basePath,'UPDOWNstates');
    return
end

%% Collect all the Necessary Pieces of Information

%Sleep Scoring (for NREM)
[SleepState] = bz_LoadStates(basePath,'SleepState');
if isempty(SleepState) && isempty(NREMInts)
    button = questdlg(['SleepState.states.mat does not exist, '...
        'would you like to run SleepScoreMaster?'],...
        'DetectLFPeakDOWN Needs NREM',...
        'Yes','No, use all timepoints','Cancel','Yes');
    switch button
        case 'Yes'
            SleepState = SleepScoreMaster(basePath);
            display('Please double check quality of sleep scoring in the StateScoreFigures folder')
            NREMInts = SleepState.ints.NREMstate;
        case 'No, use all timepoints'
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
    SWChann = input('Which channel shows the most robust (positive polarity) slow waves?');
elseif isempty(SWChann)
    SWChann = SleepState.detectorparams.SWchannum;
end
lfp = bz_GetLFP(SWChann,'basepath',basePath);

%Spikes in the CTX Spike Groups
spikes = bz_GetSpikes('basepath',basePath,'spikeGroups',CTXSpikeGroups);


%% Filter the LFP
display('Filtering LFP')
filterbounds = [0.5 6]; %heuristically defined.  room for improvement here.
deltaLFP = bz_Filter(lfp,'passband',filterbounds,'filter','fir1','order',1);

%% Find peaks in the filtered LFP
normLFP = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

[peakheights,LFPeaks] = findpeaks(normLFP,deltaLFP.timestamps,'MinPeakHeight',lowerpeakthresh,'MinPeakDistance',peakdist);
[LFPeaks,keepPeaks] = RestrictInts(LFPeaks,NREMInts);LFPeaks=LFPeaks(:,1);
LFPeakheight = peakheights(keepPeaks);

%% Check cell number and put all the spikes together
numcells = length(spikes.times);
if numcells < goodcellnumber
    warning(['You only have ',num2str(numcells),' in this recording... >',...
        num2str(goodcellnumber),' are recommended for good SW detection'])
end
allspikes = sort(cat(1,spikes.times{:}));


%% Peak threshold for putative SWs

putSWs = LFPeaks(LFPeakheight>=peakthresh);
putSWPeakHeights = LFPeakheight(LFPeakheight>=peakthresh);

%% Calculate Binned Population Rate
display('Binning Spikes')
dt = 0.005; %dt = 5ms
overlap = 6; %Overlap = 6 dt
winsize = dt*overlap; %meaning windows are 30ms big
[spikemat,t_spkmat,spindices] = SpktToSpkmat(spikes.times, [], dt,overlap);
synchmat = sum(spikemat>0,2);
ratemat = sum(spikemat,2);

%% NREM spike rate histogram

[tspike_NREM,tidx_NREM] = RestrictInts(t_spkmat,NREMInts);
tspike_NREM = tspike_NREM(:,1);
synchmat_NREM = synchmat(tidx_NREM);
ratemat_NREM = ratemat(tidx_NREM);

%% Find OFF states

minoffspikes = 1;
OFFidx = synchmat<=minoffspikes;
OFFints = IDXtoINT(OFFidx);
OFFints = t_spkmat(OFFints{1});
OFFints = RestrictInts(OFFints,NREMInts);


%% Find DOWN states: Delta peaks in OFF ints
%maxinterdown = 15;
minDOWNdur = 0.03;
OFFints(OFFints(:,2)-OFFints(:,1)<minDOWNdur,:) = [];
[~,keepdelta,keepOFF] = RestrictInts(putSWs,OFFints);
DOWNints = OFFints(keepOFF,:);
DOWNpeaks = putSWs(keepdelta);


DOWNpeakmag = putSWPeakHeights(keepdelta);
DOWNdur = diff(DOWNints,[],2);

%Remove DOWNs that are shorter than the threshold
removepeaks = DOWNdur<minOFF;
DOWNints(removepeaks,:) = [];DOWNdur(removepeaks) = [];
DOWNpeaks(removepeaks) = []; DOWNpeakmag(removepeaks) = [];

%Calculate UPs and remove those that aren't in detectionints
UPints = [DOWNints(1:end-1,2) DOWNints(2:end,1)];
UPints = RestrictInts(UPints,NREMInts);

UPdur = diff(UPints,[],2);



%%    
if SHOWFIG  
%% For Figure
display('Calculating PETH by Peak for Figure')
win = 1;
reltime = [];
spkpeakheight = [];
for pp = 1:length(LFPeaks)
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
ratemat_bypeakmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);

    
%% CCG of SLow Waves
% CCGvec = [DOWNpeaks,ones(size(DOWNpeaks));...
%     DOWNints(:,1),2.*ones(size(DOWNints(:,1)));...
%     DOWNints(:,2),3.*ones(size(DOWNints(:,2)))];
%With other slow waves
[SW_CCG,t_CCG] = CCG(DOWNpeaks,ones(size(DOWNpeaks)),'binSize',0.02);

%With Spikes
CCGvec = [DOWNpeaks,ones(size(DOWNpeaks));...
    allspikes,2.*ones(size(allspikes))];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02);
%SWspike_CCG = SWspike_CCG(:,1,2);

%%
winsize = 10; %s
samplewin = randsample(LFPeaks,1)+winsize.*[-0.5 0.5];

ratecolor = makeColorMap([1 1 1],[0.8 0 0],[0 0 0]);
figure

subplot(8,2,2)
    hist(ratemat_NREM,0:numcells)
    hold on
    plot([0 0],get(gca,'ylim'),'k-')
    xlim([0 20])
subplot(8,2,4)
    hist(synchmat_NREM,0:numcells)
    hold on
    plot([0 0],get(gca,'ylim'),'k-')
    xlim([0 20])

    
subplot(4,2,3)
    imagesc(timebins,magbins,ratemat_bypeakmag')
    hold on
    plot(get(gca,'xlim'),peakthresh.*[1 1],'k--')
    plot([0 0],get(gca,'ylim'),'k-')
    colormap(ratecolor)
    colorbar('location','east')
    xlim([-0.8 0.8])
    axis xy
    xlabel('t (relative to delta peak)');ylabel({'Delta Peak Amplitude', '(modZ)'})
    
subplot(8,2,1)
    bar(t_CCG,SW_CCG,'facecolor','g','FaceAlpha',0.2)
    hold on
    plot([0 0],get(gca,'ylim'),'k-')
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Delta Peaks')
    
subplot(8,2,3)
    bar(t_CCG,SWspike_CCG(:,1,2)./length(DOWNpeaks),'facecolor',[0.5 0.5 0.5])
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Spike Rate');
    
subplot(8,2,6)
    hist(log10(UPdur),40)
    LogScale('x',10)
subplot(8,2,8)
    hist(log10(DOWNdur),40)
    LogScale('x',10)
   
    
subplot(6,1,4)
    plot(lfp.timestamps,lfp.data,'k')
    hold on
    box off
    plot(lfp.timestamps(ismember(lfp.timestamps,LFPeaks)),lfp.data(ismember(lfp.timestamps,LFPeaks)),'r.')
    plot(lfp.timestamps(ismember(lfp.timestamps,DOWNpeaks)),lfp.data(ismember(lfp.timestamps,DOWNpeaks)),'g.')
    
        %OFF patches
        y = get(gca,'Ylim');
        patch([OFFints(:,1) OFFints(:,2) OFFints(:,2) OFFints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(OFFints(:,1)),1)',...
            'r','FaceAlpha',0.1,'EdgeColor','none');
    
        %DOWN patches
        y = get(gca,'Ylim');
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
    
    xlim(samplewin)
    set(gca,'xticklabel',[]);  set(gca,'ytick',[])
    ylabel('LFP')
    
subplot(3,1,3)
    bar(tspike_NREM,synchmat_NREM,'facecolor',[0.5 0.5 0.5])
    hold on
    ylim([-8 length(unique(spikes.spindices(:,2)))+1])
        %OFF patches
        y = get(gca,'Ylim');
        patch([OFFints(:,1) OFFints(:,2) OFFints(:,2) OFFints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(OFFints(:,1)),1)',...
            'r','FaceAlpha',0.1,'EdgeColor','none');    
        %DOWN patches
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
     
        
    plot(deltaLFP.timestamps,normLFP-peakthresh,'b','Linewidth',1)    
    plot(spikes.spindices(:,1),spikes.spindices(:,2),'k.','MarkerSize',4)
    box off    
    xlim(samplewin)
    ylabel('Cells');xlabel('t (s)')
    legend('Spike Synchrony','OFF States','DOWN States','Filtered LFP')
    
NiceSave('SlowOscillation',figfolder,baseName)

end


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
UPDOWNstates.detectiondate = today('datetime');

if SAVEMAT
    save(savefile,'UPDOWNstates')
end

end

