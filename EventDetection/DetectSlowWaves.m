function [ SlowWaves ] = DetectSlowWaves( basePath,varargin)
%[UPDOWNstates] = DetectSlowWave(basePath) detects neocortical slow
%waves using a combination of a positive deflection in the LFP (delta wave)
%and a dip in gamma power.
%
%INPUTS
%   basePath  
%   (input options not yet implemented... all set to default)
%   'NREMInts'  -Interval of times for NREM 
%               -(Default: loaded from SleepState.states.mat)
%   'SWChan'   -Channel with the most robust (positively deflecting) Slow
%               Waves. (0-Indexing a la neuroscope). 
%               can try 'autoselect'
%               'useold' to use channel from existing SlowWave.states.mat
%               (Default: use SWChannel from sleep scoring)
%   'CTXSpkGroups' -Spike groups that are in the cortex...  default: all
%   'CTXChans' -LFP channels that are in the cortex...  default: all
%   'detectionparms' -Structure of detection parameters (needs update)
%       .minOFF         (default: 0.025s)
%       .mininterOFF    (default: 0.03s)
%       .peakdist       (default: 0.05s)
%       .peakthresh     (default: 2.5 stds - modified Z score)
%   'SHOWFIG'   -true/false show a quality control figure (default: true)
%   'saveMat'   -logical (default=true) to save in buzcode format
%   'forceReload' -logical (default: false) to redetect (add option to use
%                   old parameters like channels...)
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
p = inputParser;
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'showFig',true,@islogical);
addParameter(p,'SWChan',[]);
addParameter(p,'NREMInts',[]);
addParameter(p,'CTXChans','all');
parse(p,varargin{:})

FORCEREDETECT = p.Results.forceReload;
SAVEMAT = p.Results.saveMat;
SHOWFIG = p.Results.showFig;
SWChann = p.Results.SWChan;
NREMInts = p.Results.NREMInts;
CTXChans = p.Results.CTXChans;


%Defaults
if ~exist('basePath','var')
    basePath = pwd;
end

%Parms
% DELTApeakthresh = 2.2;
% DELTAwinthresh = 1;
% GAMMAdipthresh = 1.25;
% GAMMAwinthresh = 1;

minwindur = 0.04;
joinwindur = 0.005;

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.SlowWaves.events.mat']);

if exist(savefile,'file') && ~FORCEREDETECT
    display(['Slow Oscillation already Detected, loading ',baseName,'.SlowWave.states.mat'])
    SlowWaves = bz_LoadStates(basePath,'SlowWave');
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
    CHANSELECT = 'userinput';
    SWChann = input('Which channel shows the most robust (positive polarity) slow waves?');
elseif isempty(SWChann) %&& ismember(SleepState.detectorparams.SWchannum,CTXChans)
    CHANSELECT = 'SleepScoreSWchan';
    SWChann = SleepState.detectorparams.SWchannum;
elseif strcmp(SWChann,'manualselect')
    CHANSELECT = 'manual';
    %run SW channel selection routine: subfunction below (ManChanSelect)
elseif strcmp(SWChann,'autoselect')
    CHANSELECT = 'auto';
    %run SW channel selection routine: subfunction below (AutoChanSelect)
    SWChann = AutoChanSelect(CTXChans,basePath,NREMInts);
elseif strcmp(SWChann,'useold')
    SlowWaves = bz_LoadStates(basePath,'SlowWave');
    CHANSELECT = SlowWaves.detectorinfo.detectionparms.CHANSELECT;
    SWChann = SlowWaves.detectorinfo.detectionparms.SWchannel;
else
    CHANSELECT = 'userinput';
end
lfp = bz_GetLFP(SWChann,'basepath',basePath);

%Spikes in the CTX Spike Groups
spikes = bz_GetSpikes('basepath',basePath,'region','CTX');
allspikes = sort(cat(1,spikes.times{:}));
%assumes region 'CTX'.... update this maybe putting in local region 
%spikegroups is the better way to go

%% Filter the LFP: delta
display('Filtering LFP')
deltafilterbounds = [0.5 6]; %heuristically defined.  room for improvement here.
deltaLFP = bz_Filter(lfp,'passband',deltafilterbounds,'filter','fir1','order',1);
deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

%% Filter and get power of the LFP: gamma
gammafilter = [100 inf]; %high pass >80Hz (previously (>100Hz)
gammasmoothwin = 0.08; %window for smoothing gamma power 
gammaLFP = bz_Filter(lfp,'passband',gammafilter,'filter','fir1','order',4);
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(gammasmoothwin.*gammaLFP.samplingRate),'moving' );
gammaLFP.normamp = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,gammaLFP.samplingRate);

%% Determine Thresholds and find windows around delta peaks/gamma dips
[thresholds,threshfigs] = DetermineThresholds(deltaLFP,gammaLFP,spikes,NREMInts);

%Find peaks in the delta LFP
[DELTApeaks,DELTAwins,DELTApeakheight] = FindPeakInWin(deltaLFP.normamp,deltaLFP.timestamps,...
    thresholds.DELTApeakthresh,thresholds.DELTAwinthresh,minwindur,joinwindur);
[DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);
DELTApeakheight = DELTApeakheight(keepPeaks);  DELTAwins = DELTAwins(keepPeaks,:);

%Find dips in the gamma power
[GAMMAdips,GAMMAwins,GAMMAdipdepth] = FindPeakInWin(-gammaLFP.normamp,gammaLFP.timestamps,...
    thresholds.GAMMAdipthresh,thresholds.GAMMAwinthresh,minwindur,joinwindur);
[GAMMAdips,keepPeaks] = RestrictInts(GAMMAdips,NREMInts);
GAMMAdipdepth = GAMMAdipdepth(keepPeaks);  GAMMAwins = GAMMAwins(keepPeaks,:);

%% Merge gamma/delta windows to get Slow Waves, UP/DOWN states
[ DOWNints,mergedidx ] = MergeSeparatedInts( [DELTAwins;GAMMAwins]);

if any(DOWNints(2:end,1)-DOWNints(1:end-1,2)<=0)
    display('Merge Error... why?')
end

%Keep only those windows in which DELTA/GAMMA were merged togehter (note
%this only works if joining happened previously... otherwise could keep
%windows where two delta and/or two gamma were joined; FindPeakInWin joins
numwins = cellfun(@length,mergedidx);
DOWNints = DOWNints(numwins>=2,:);
%Get the SW peak magnitude
mergedidx = mergedidx(numwins>=2); %keep only the indices from those that were merged into SWs
mergeddeltaidx = cellfun(@(X) X(X<=length(DELTAwins)),mergedidx,'UniformOutput',false); %keep delta
[SWpeakmag,peakidx] = cellfun(@(X) max(DELTApeakheight(X)),mergeddeltaidx,'UniformOutput',false); %Pick the larger of the peaks in each SW
SWpeaks = cellfun(@(X,Y) DELTApeaks(X(Y)),mergeddeltaidx,peakidx);
SWpeakmag = [SWpeakmag{:}]';

%UP and DOWN
DOWNdur = diff(DOWNints,[],2);
%DL: duration threshold already accounted for
% %Remove DOWNs that are shorter than the threshold
% removepeaks = DOWNdur<minOFF;
% DOWNints(removepeaks,:) = [];DOWNdur(removepeaks) = [];
% DOWNpeaks(removepeaks) = []; DOWNpeakmag(removepeaks) = [];

%Calculate UPs and remove those that aren't in detectionints
UPints = [DOWNints(1:end-1,2) DOWNints(2:end,1)];
UPints = RestrictInts(UPints,NREMInts);
UPdur = diff(UPints,[],2);



%%
winsize = 7; %s
samplewin = randsample(DELTApeaks,1)+winsize.*[-0.5 0.5];
sampleIDX = lfp.timestamps>=samplewin(1) & lfp.timestamps<=samplewin(2);

figure
subplot(2,1,1)
plot(deltaLFP.timestamps(sampleIDX),gammaLFP.normamp(sampleIDX),'g')
hold on
plot(deltaLFP.timestamps(sampleIDX),deltaLFP.normamp(sampleIDX),'r')
plot(GAMMAdips,-GAMMAdipdepth,'go')
plot(DELTApeaks,DELTApeakheight,'ro')
plot(SWpeaks,SWpeakmag,'ko')
plot(GAMMAwins',-thresholds.GAMMAwinthresh.*ones(size(GAMMAwins')),'g')
plot(DELTAwins',thresholds.DELTAwinthresh.*ones(size(DELTAwins')),'r')
xlim(samplewin)
subplot(2,1,2)
plot(lfp.timestamps(sampleIDX),lfp.data(sampleIDX),'k')
hold on
plot(DOWNints',zeros(size(DOWNints')),'k')
xlim(samplewin)
%plot(gammaLFP.normamp,deltaLFP.normamp,'k.')
xlabel('Gamma Power');ylabel('Delta LFP')





%%   Stuff for the detection figure
if SHOWFIG  

%% Calculate Binned Population Rate
display('Binning Spikes')
dt = 0.005; %dt = 5ms
overlap = 8; %Overlap = 8 dt
winsize = dt*overlap; %meaning windows are 40ms big (previously 30)
[spikemat,t_spkmat,spindices] = SpktToSpkmat(spikes.times, [], dt,overlap);
synchmat = sum(spikemat>0,2);
ratemat = sum(spikemat,2);

%NREM spike rate histogram
[tspike_NREM,tidx_NREM] = RestrictInts(t_spkmat,NREMInts);
tspike_NREM = tspike_NREM(:,1);
synchmat_NREM = synchmat(tidx_NREM);
ratemat_NREM = ratemat(tidx_NREM);

%% CCG of SLow Waves and spikes
[SW_CCG,t_CCG] = CCG(SWpeaks,ones(size(SWpeaks)),'binSize',0.02);

%With Spikes
CCGvec = [SWpeaks,ones(size(SWpeaks));...
    allspikes,2.*ones(size(allspikes))];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02);

%%
winsize = 10; %s
samplewin =randsample(SWpeaks,1)+winsize.*[-0.5 0.5];
sampleIDX = lfp.timestamps>=samplewin(1) & lfp.timestamps<=samplewin(2);

ratecolor = makeColorMap([1 1 1],[0.8 0 0],[0 0 0]);

figure

    gammafig = copyobj(threshfigs.GAMMArate,gcf);
    subplot(4,2,4,gammafig)
    
    deltafig = copyobj(threshfigs.DELTArate,gcf);
    subplot(4,2,3,deltafig)

    
subplot(8,2,1)
    bar(t_CCG,SW_CCG,'facecolor','g','FaceAlpha',0.2)
    hold on
    plot([0 0],get(gca,'ylim'),'k-')
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Delta Peaks')
    
subplot(8,2,3)
    bar(t_CCG,SWspike_CCG(:,1,2)./length(SWpeaks),'facecolor',[0.5 0.5 0.5])
    xlim([-0.8 0.8])
    set(gca,'xticklabel',[]);
    ylabel('Spike Rate');
    
subplot(8,2,2)
    hist(log10(UPdur),40)
    LogScale('x',10)
subplot(8,2,4)
    hist(log10(DOWNdur),40)
    LogScale('x',10)
   
    
subplot(6,1,4)
    plot(lfp.timestamps(sampleIDX),lfp.data(sampleIDX),'k')
    axis tight
    hold on
    box off
    plot(lfp.timestamps(ismember(lfp.timestamps,SWpeaks)),lfp.data(ismember(lfp.timestamps,SWpeaks)),'r.')
    plot(lfp.timestamps(ismember(lfp.timestamps,SWpeaks)),lfp.data(ismember(lfp.timestamps,SWpeaks)),'g.')
    
        %DOWN patches
        y = get(gca,'Ylim');
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
    
    xlim(samplewin)
    set(gca,'xticklabel',[]);  set(gca,'ytick',[])
    ylabel({'LFP',['Chan: ',num2str(SWChann)]})
    
subplot(3,1,3)
    bar(tspike_NREM,synchmat_NREM,'facecolor',[0.5 0.5 0.5])
    hold on
    ylim([-8 length(unique(spikes.spindices(:,2)))+1])
        %DOWN patches
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
     
   % plot(deltaLFP.timestamps,deltaLFP.normamp-peakthresh,'b','Linewidth',1)    
    plot(deltaLFP.timestamps(sampleIDX),gammaLFP.normamp(sampleIDX)+thresholds.GAMMAwinthresh,'g','Linewidth',1) 
    plot(deltaLFP.timestamps(sampleIDX),deltaLFP.normamp(sampleIDX)-thresholds.DELTAwinthresh,'b','Linewidth',1) 
    plot(spikes.spindices(:,1),spikes.spindices(:,2),'k.','MarkerSize',4)
    box off    
    xlim(samplewin)
    ylabel('Cells');xlabel('t (s)')
    legend('Spike Synchrony','OFF States','DOWN States','Filtered LFP')
    
NiceSave('SlowOscillation',figfolder,baseName)

end


%% Ouput in .event.mat format
%Needs to be updated
detectionparms.SWchannel = SWChann;
detectionparms.CHANSELECT = CHANSELECT;
detectionparms.CTXChans = CTXChans;
detectionparms.thresholds = thresholds;

SlowWaves.ints.UP = UPints;
SlowWaves.ints.DOWN = DOWNints;
SlowWaves.timestamps = SWpeaks;
SlowWaves.SWpeakmag = SWpeakmag;
SlowWaves.detectorinfo.detectorname = 'DetectSlowWaves';
SlowWaves.detectorinfo.detectionparms = detectionparms;
SlowWaves.detectorinfo.detectiondate = today('datetime');

if SAVEMAT
    save(savefile,'SlowWave')
end

display('Slow Wave Detection: COMPLETE!')
end








%% Channel Selection Functions
function usechan = AutoChanSelect(trychans,basePath,NREMInts)
    display('Detecting best channel for slow wave detection...')
    baseName = bz_BasenameFromBasepath(basePath);
    figfolder = fullfile(basePath,'DetectionFigures');
    %Exclude badchannels
    par = bz_getSessionInfo(basePath);
    if strcmp(trychans,'all') 
        trychans = [par.SpkGrps(:).Channels];
    end
    if isfield(par,'badchannels')
    	trychans = setdiff(trychans,par.badchannels);
    end
    
    for cc = 1:length(trychans)
        display(['Trying Channel ',num2str(cc),' of ',num2str(length(trychans))])
        %Load the LFPs
        chanlfp = bz_GetLFP(trychans(cc),'basepath',basePath);
        %Filter in gamma
        gammafilter = [100 inf];
        trygammaLFP = bz_Filter(chanlfp,'passband',gammafilter,'filter','fir1','order',4);

        %Restrict to NREM only - could also use intervals above to do this....
        [chanlfp.timestamps,inNREMidx] = RestrictInts(chanlfp.timestamps,NREMInts);
        chanlfp.data = chanlfp.data(inNREMidx,:);
        trygammaLFP.amp = trygammaLFP.amp(inNREMidx,:);
        %Best channel is the one in which gamma is most anticorrelated with the
        %LFP - i.e. DOWN states (positive LFP) have low gamma power
    
        gammaLFPcorr(cc) = corr(single(chanlfp.data),trygammaLFP.amp,'type','spearman');
        
        %Save a small window for plotting
        if ~exist('samplewin','var')
            winsize = 4; %s
            samplewin =randsample(chanlfp.timestamps,1)+winsize.*[-0.5 0.5];
            sampleIDX = chanlfp.timestamps>=samplewin(1) & chanlfp.timestamps<=samplewin(2);
            alllfp.timestamps = chanlfp.timestamps(sampleIDX);
        end
        alllfp.data(:,cc) = chanlfp.data(sampleIDX);
        alllfp.channels(cc) = chanlfp.channels;
    end
    
    [~,usechanIDX] = min(gammaLFPcorr);
    usechan = trychans(usechanIDX);
    
    display(['Selected Channel: ',num2str(usechan)])
    
    [~,sortcorr] = sort(gammaLFPcorr);
    %% Figure

    
    figure
    subplot(2,2,1)
        hist(gammaLFPcorr)
        xlabel('Gamma-LFP Correlation')
%     subplot(2,2,3)
%         plot(trygammaLFP.amp(:,usechanIDX),chanlfp.data(:,usechanIDX),'k.')
%         xlabel('Gamma Power');ylabel('Raw LFP')
    subplot(6,2,4:2:12)
        bz_MultiLFPPlot(alllfp,'channels',trychans(sortcorr),'timewin',samplewin)
        ylabel('Slow Wave Channel (Worse<--------->Better)')
        
    subplot(6,2,2)
        plot(alllfp.timestamps,alllfp.data(:,usechanIDX),'k')
        axis tight
        box off
        
    NiceSave('SlowWaveChannelSelect',figfolder,baseName)
end



%% Threshold Determination Function
function [thresholds,threshfigs] = DetermineThresholds(deltaLFP,gammaLFP,spikes,NREMInts)
    display('Determining delta/gamma thresholds for detection...')
    minwindur = 0.04; %should pass through... but not super important
    
    %Find peaks in delta, gamma power
    [peakheights,DELTApeaks] = findpeaks(deltaLFP.normamp,deltaLFP.timestamps,'MinPeakHeight',0.25,'MinPeakDistance',minwindur);
    [DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);
    DELTAPeakheight = peakheights(keepPeaks);
    [~,DELTApeakIDX] = ismember(DELTApeaks,deltaLFP.timestamps);

    [peakheights,GAMMAdips] = findpeaks(-gammaLFP.normamp,gammaLFP.timestamps,'MinPeakHeight',0.6,'MinPeakDistance',minwindur);
    [GAMMAdips,keepPeaks] = RestrictInts(GAMMAdips,NREMInts);
    GAMMAdipdepth = peakheights(keepPeaks);
    [~,GAMMAdipIDX] = ismember(GAMMAdips,gammaLFP.timestamps);
    
    %% Get the spike PETH around delta/gamma peaks
    display('Calculating PETH by Peak For Threshold Calibration')
    win = 1;
    numchecks = 15000;

    numtimebins = 200;
    nummagbins = 30;
    timebins = linspace(-win,win,numtimebins);
    
    numcells = length(spikes.times);
    allspikes = sort(cat(1,spikes.times{:}));
    
    %DELTA
    display('DELTA...')
    reltime = [];
    spkpeakheight = [];
    if length(DELTApeaks)>numchecks
        sampleDELTA = randsample(length(DELTApeaks),numchecks);
    else
        sampleDELTA = 1:length(DELTApeaks);
    end
    nearpeakdelta = zeros(2.*win.*deltaLFP.samplingRate+1,nummagbins);
    DELTAmagbins = linspace(min(DELTAPeakheight),min([max(DELTAPeakheight)-1,8]), nummagbins);
    for pp = 1:length(sampleDELTA)
        ss = sampleDELTA(pp);

        %Spikes around the delta
        nearpeakspikes = allspikes >= DELTApeaks(ss)-win & allspikes <= DELTApeaks(ss)+win;
        reltime = [reltime; allspikes(nearpeakspikes)-DELTApeaks(ss)];
        spkpeakheight = [spkpeakheight; DELTAPeakheight(ss).*ones(sum(nearpeakspikes),1)];

        %Delta around the delta
        [~,groupidx] = min(abs(DELTAPeakheight(ss)-DELTAmagbins));
        nearpeakdelta(:,groupidx) =nearpeakdelta(:,groupidx) + ...
            deltaLFP.normamp(DELTApeakIDX(ss)+[-1.*deltaLFP.samplingRate:deltaLFP.samplingRate]);
    end

    peakmagdist = hist(DELTAPeakheight(sampleDELTA),DELTAmagbins);
    spikehitmat = hist3([reltime,spkpeakheight],{timebins,DELTAmagbins});
    ratemat_byDELTAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);
    deltapower_byDELTAmag = bsxfun(@(A,B) A./B,nearpeakdelta,peakmagdist);
    
    %GAMMA
    display('GAMMA...')
    reltime = [];
    spkpeakheight = [];
    if length(GAMMAdips)>numchecks
        sampleGAMMA = randsample(length(GAMMAdips),numchecks); %Should really try to get even# in each bin...
    else
        sampleGAMMA = 1:length(GAMMAdips);
    end
    GAMMAmagbins = linspace(min(GAMMAdipdepth),min([max(GAMMAdipdepth)-0.5,8]),nummagbins);
    neardipgamma = zeros(2.*win.*gammaLFP.samplingRate+1,nummagbins);
    for pp = 1:length(sampleGAMMA)
        ss = sampleGAMMA(pp);
        nearpeakspikes = allspikes >= GAMMAdips(ss)-win & allspikes <= GAMMAdips(ss)+win;
        reltime = [reltime; allspikes(nearpeakspikes)-GAMMAdips(ss)];
        spkpeakheight = [spkpeakheight; GAMMAdipdepth(ss).*ones(sum(nearpeakspikes),1)];

         %Gamma around the gamma
        [~,groupidx] = min(abs(GAMMAdipdepth(ss)-GAMMAmagbins));
        neardipgamma(:,groupidx) =neardipgamma(:,groupidx) + ...
            gammaLFP.normamp(GAMMAdipIDX(ss)+[-1.*gammaLFP.samplingRate:gammaLFP.samplingRate]);
    end

    peakmagdist = hist(GAMMAdipdepth(sampleGAMMA),GAMMAmagbins);
    spikehitmat = hist3([reltime,spkpeakheight],{timebins,GAMMAmagbins});
    ratemat_byGAMMAmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);
    gammapower_byGAMMAmag = bsxfun(@(A,B) A./B,neardipgamma,peakmagdist);
    %% Find Thresholds as values where spike rate falls below threshold
    ratethresh = 0.2;

    DELTAbox=bwmorph(ratemat_byDELTAmag<ratethresh,'close');
    DELTAbox=bwmorph(DELTAbox,'open');
    DELTAbox=bwboundaries(DELTAbox); DELTAbox = DELTAbox{1};
    DELTAbox(DELTAbox(:,2)==max(DELTAbox(:,2)),:)=[];
    DELTApeakthresh = DELTAmagbins(min(DELTAbox(:,2)));
    scalefactor = length(nearpeakdelta)./numtimebins; %convert number of bins for LFP and spikes
    DELTAwinthresh = deltapower_byDELTAmag(sub2ind(size(deltapower_byDELTAmag),...
        round(DELTAbox(:,1).*scalefactor),DELTAbox(:,2)));
    DELTAwinthresh = mean(DELTAwinthresh);

    GAMMAbox=bwmorph(ratemat_byGAMMAmag<ratethresh,'close');
    GAMMAbox=bwmorph(GAMMAbox,'open');
    if sum(GAMMAbox(:))== 0
        display('No DOWN around gamma dip.... perhaps adjust rate threshold or pick another channel?')
        GAMMAwinthresh = 1.2;
        GAMMAwinthresh = 1;
        GAMMAbox = [1 1];
    else
    GAMMAbox=bwboundaries(GAMMAbox); GAMMAbox = GAMMAbox{1};
    GAMMAbox(GAMMAbox(:,2)==max(GAMMAbox(:,2)),:)=[];
    GAMMAdipthresh = GAMMAmagbins(min(GAMMAbox(:,2)));
    scalefactor = length(neardipgamma)./numtimebins; %convert number of bins for LFP and spikes
    GAMMAwinthresh = gammapower_byGAMMAmag(sub2ind(size(gammapower_byGAMMAmag),...
        round(GAMMAbox(:,1).*scalefactor),GAMMAbox(:,2)));
    GAMMAwinthresh = -mean(GAMMAwinthresh);
    end
    
    thresholds.DELTApeakthresh = DELTApeakthresh;
    thresholds.DELTAwinthresh = DELTAwinthresh;
    thresholds.GAMMAdipthresh = GAMMAdipthresh;
    thresholds.GAMMAwinthresh = GAMMAwinthresh;
    thresholds.ratethresh = ratethresh;

    %%
    figure
     threshfigs.DELTArate = subplot(2,2,3);
        imagesc(timebins,DELTAmagbins,ratemat_byDELTAmag')
        hold on
        plot(timebins(DELTAbox(:,1)),DELTAmagbins(DELTAbox(:,2)),'r.')
        plot(get(gca,'xlim'),DELTApeakthresh.*[1 1],'k--')
        xlabel('t (relative to SW peak)');ylabel({'SW Peak Amplitude', '(modZ)'})
        axis xy

    threshfigs.GAMMArate = subplot(2,2,4);
        imagesc(timebins,GAMMAmagbins,ratemat_byGAMMAmag')
        hold on
        plot(timebins(GAMMAbox(:,1)),GAMMAmagbins(GAMMAbox(:,2)),'r.')
        plot(get(gca,'xlim'),GAMMAdipthresh.*[1 1],'k--')
        xlabel('t (relative to GA Dip)');ylabel({'GA Dip Amplitude', '(modZ)'})

    threshfigs.DELTApower = subplot(2,2,1);
        imagesc(timebins,DELTAmagbins,deltapower_byDELTAmag')
        hold on
        plot(timebins(DELTAbox(:,1)),DELTAmagbins(DELTAbox(:,2)),'r.')
        plot(get(gca,'xlim'),DELTApeakthresh.*[1 1],'k--')
        colorbar
        xlabel('t (relative to SW peak)');ylabel({'SW Peak Amplitude', '(modZ)'})
        axis xy

    threshfigs.GAMMApower = subplot(2,2,2);
        imagesc(timebins,GAMMAmagbins,gammapower_byGAMMAmag')
        hold on
        plot(timebins(GAMMAbox(:,1)),GAMMAmagbins(GAMMAbox(:,2)),'r.')
        plot(get(gca,'xlim'),GAMMAdipthresh.*[1 1],'k--')
        xlabel('t (relative to GA Dip)');ylabel({'GA Dip Amplitude', '(modZ)'})
        colorbar

end
    