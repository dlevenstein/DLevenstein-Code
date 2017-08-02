function [ SlowWave ] = DetectSlowWaves( basePath,varargin)
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
DELTApeakthresh = 2.2;
DELTAwinthresh = 1;
GAMMAdipthresh = 1.25;
GAMMAwinthresh = 1;

minwindur = 0.04;
joinwindur = 0.005;

%% File Management
baseName = bz_BasenameFromBasepath(basePath);
figfolder = fullfile(basePath,'DetectionFigures');
savefile = fullfile(basePath,[baseName,'.SlowWave.states.mat']);

if exist(savefile,'file') && ~FORCEREDETECT
    display(['Slow Oscillation already Detected, loading ',baseName,'.SlowWave.states.mat'])
    SlowWave = bz_LoadStates(basePath,'SlowWave');
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
else
    CHANSELECT = 'userinput';
end
lfp = bz_GetLFP(SWChann,'basepath',basePath);

%Spikes in the CTX Spike Groups
spikes = bz_GetSpikes('basepath',basePath,'region','CTX');
%assumes region 'CTX'.... update this maybe putting in local region 
%spikegroups is the better way to go

%% Check cell number and put all the spikes together
 numcells = length(spikes.times);
% if numcells < goodcellnumber
%     warning(['You only have ',num2str(numcells),' in this recording... >',...
%         num2str(goodcellnumber),' are recommended for good SW detection'])
% end
allspikes = sort(cat(1,spikes.times{:}));

%% Filter the LFP: delta
display('Filtering LFP')
deltafilterbounds = [0.5 6]; %heuristically defined.  room for improvement here.
deltaLFP = bz_Filter(lfp,'passband',deltafilterbounds,'filter','fir1','order',1);

%% Find peaks in the delta LFP
normDELTA = NormToInt(deltaLFP.data,'modZ',NREMInts,deltaLFP.samplingRate);

% [peakheights,DELTApeaks] = findpeaks(normDELTA,deltaLFP.timestamps,'MinPeakHeight',lowerpeakthresh,'MinPeakDistance',peakdist);
% [DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);DELTApeaks=DELTApeaks(:,1);
% DELTAPeakheight = peakheights(keepPeaks);

[DELTApeaks,DELTAwins,DELTApeakheight] = FindPeakInWin(normDELTA,deltaLFP.timestamps,DELTApeakthresh,DELTAwinthresh,minwindur,joinwindur);
[DELTApeaks,keepPeaks] = RestrictInts(DELTApeaks,NREMInts);
DELTApeakheight = DELTApeakheight(keepPeaks);  DELTAwins = DELTAwins(keepPeaks,:);

%Peak threshold for putative SWs
% putSWs = DELTApeaks(DELTAPeakheight>=peakthresh);
% putSWPeakHeights = DELTAPeakheight(DELTAPeakheight>=peakthresh);

%% Filter and get power of the LFP: gamma
gammafilter = [100 inf]; %high pass >80Hz (previously (>100Hz)
gammasmoothwin = 0.07; %window for smoothing gamma power 
gammaLFP = bz_Filter(lfp,'passband',gammafilter,'filter','fir1','order',4);
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(gammasmoothwin.*gammaLFP.samplingRate),'moving' );

%% Find dips in the gamma power
normGAMMA = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,gammaLFP.samplingRate);

[GAMMAdips,GAMMAwins,GAMMAdipdepth] = FindPeakInWin(-normGAMMA,gammaLFP.timestamps,GAMMAdipthresh,GAMMAwinthresh,minwindur,joinwindur);
[GAMMAdips,keepPeaks] = RestrictInts(GAMMAdips,NREMInts);
GAMMAdipdepth = GAMMAdipdepth(keepPeaks);  GAMMAwins = GAMMAwins(keepPeaks,:);
%
%use width in findpeaks to get window below threshold? same for delta?

%% Merge gamma/delta windows
[ DOWNints,mergedidx ] = MergeSeparatedInts( [DELTAwins;GAMMAwins]);

if any(DOWNints(2:end,1)-DOWNints(1:end-1,2)<=0)
    display('Merge Error... why?')
end

%Keep only those windows in which DELTA/GAMMA were merged togehter (note
%this only works if joining happened previously... otherwise could keep
%windows where two delta and/or two gamma were joined;
numwins = cellfun(@length,mergedidx);
DOWNints = DOWNints(numwins>=2,:);
%Get the SW peak magnitude
mergedidx = mergedidx(numwins>=2); %keep only the indices from those that were merged into SWs
mergeddeltaidx = cellfun(@(X) X(X<=length(DELTAwins)),mergedidx,'UniformOutput',false); %keep delta
[SWpeakmag,peakidx] = cellfun(@(X) max(DELTApeakheight(X)),mergeddeltaidx,'UniformOutput',false); %Pick the larger of the peaks in each SW
SWpeaks = cellfun(@(X,Y) DELTApeaks(X(Y)),mergeddeltaidx,peakidx);
SWpeakmag = [SWpeakmag{:}]';

%% UP and DOWN

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
plot(deltaLFP.timestamps(sampleIDX),normGAMMA(sampleIDX),'g')
hold on
plot(deltaLFP.timestamps(sampleIDX),normDELTA(sampleIDX),'r')
plot(GAMMAdips,-GAMMAdipdepth,'go')
plot(DELTApeaks,DELTApeakheight,'ro')
plot(SWpeaks,SWpeakmag,'ko')
plot(GAMMAwins',-GAMMAwinthresh.*ones(size(GAMMAwins')),'g')
plot(DELTAwins',DELTAwinthresh.*ones(size(DELTAwins')),'r')
xlim(samplewin)
subplot(2,1,2)
plot(lfp.timestamps(sampleIDX),lfp.data(sampleIDX),'k')
hold on
plot(DOWNints',zeros(size(DOWNints')),'k')
xlim(samplewin)
%plot(normGAMMA,normDELTA,'k.')
xlabel('Gamma Power');ylabel('Delta LFP')


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

%% Find OFF states as spiking (old)
% minoffspikes = 1;
% OFFidx = synchmat<=minoffspikes;
% OFFints = IDXtoINT(OFFidx);
% OFFints = t_spkmat(OFFints{1});
% OFFints = RestrictInts(OFFints,NREMInts);
% 
% 
% %% Find DOWN states: Delta peaks in OFF ints
% %maxinterdown = 15;
% minDOWNdur = 0.03;
% OFFints(OFFints(:,2)-OFFints(:,1)<minDOWNdur,:) = [];
% [~,keepdelta,keepOFF] = RestrictInts(putSWs,OFFints);
% DOWNints = OFFints(keepOFF,:);
% DOWNpeaks = putSWs(keepdelta);
% 
% 
% DOWNpeakmag = putSWPeakHeights(keepdelta);
% DOWNdur = diff(DOWNints,[],2);
% 
% %Remove DOWNs that are shorter than the threshold
% removepeaks = DOWNdur<minOFF;
% DOWNints(removepeaks,:) = [];DOWNdur(removepeaks) = [];
% DOWNpeaks(removepeaks) = []; DOWNpeakmag(removepeaks) = [];
% 
% %Calculate UPs and remove those that aren't in detectionints
% UPints = [DOWNints(1:end-1,2) DOWNints(2:end,1)];
% UPints = RestrictInts(UPints,NREMInts);
% 
% UPdur = diff(UPints,[],2);



%%    
if SHOWFIG  
%% For Figure

%% Spike-SW CCG histogram by peak height
% display('Calculating PETH by Peak for Figure')
% win = 1;
% reltime = [];
% spkpeakheight = [];
% for pp = 1:length(DELTApeaks)
%     nearpeakspikes = allspikes >= DELTApeaks(pp)-win & allspikes <= DELTApeaks(pp)+win;
%     reltime = [reltime; allspikes(nearpeakspikes)-DELTApeaks(pp)];
%     spkpeakheight = [spkpeakheight; DELTAPeakheight(pp).*ones(sum(nearpeakspikes),1)];
%    % keyboard
% end
% 
% 
% numtimebins = 200;
% nummagbins = 30;
% timebins = linspace(-win,win,numtimebins);
% magbins = linspace(1,8, nummagbins);
% peakmagdist = hist(DELTAPeakheight,magbins);
% 
% spikehitmat = hist3([reltime,spkpeakheight],{timebins,magbins});
% ratemat_bypeakmag = bsxfun(@(A,B) A./B./mean(diff(timebins))./numcells,spikehitmat,peakmagdist);

    
%% CCG of SLow Waves
% CCGvec = [DOWNpeaks,ones(size(DOWNpeaks));...
%     DOWNints(:,1),2.*ones(size(DOWNints(:,1)));...
%     DOWNints(:,2),3.*ones(size(DOWNints(:,2)))];
%With other slow waves
[SW_CCG,t_CCG] = CCG(SWpeaks,ones(size(SWpeaks)),'binSize',0.02);

%With Spikes
CCGvec = [SWpeaks,ones(size(SWpeaks));...
    allspikes,2.*ones(size(allspikes))];
[SWspike_CCG,t_CCG] = CCG(CCGvec(:,1),CCGvec(:,2),'binSize',0.02);
%SWspike_CCG = SWspike_CCG(:,1,2);

%%
winsize = 10; %s
samplewin =randsample(SWpeaks,1)+winsize.*[-0.5 0.5];
sampleIDX = lfp.timestamps>=samplewin(1) & lfp.timestamps<=samplewin(2);

ratecolor = makeColorMap([1 1 1],[0.8 0 0],[0 0 0]);
figure

%Add back soon
% subplot(4,2,3)
%     imagesc(timebins,magbins,ratemat_bypeakmag')
%     hold on
%     plot(get(gca,'xlim'),peakthresh.*[1 1],'k--')
%     plot([0 0],get(gca,'ylim'),'k-')
%     colormap(ratecolor)
%     colorbar('location','east')
%     xlim([-0.8 0.8])
%     axis xy
%     xlabel('t (relative to SW peak)');ylabel({'SW Peak Amplitude', '(modZ)'})
    
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
    
subplot(8,2,6)
    hist(log10(UPdur),40)
    LogScale('x',10)
subplot(8,2,8)
    hist(log10(DOWNdur),40)
    LogScale('x',10)
   
    
subplot(6,1,4)
    plot(lfp.timestamps(sampleIDX),lfp.data(sampleIDX),'k')
    axis tight
    hold on
    box off
    plot(lfp.timestamps(ismember(lfp.timestamps,SWpeaks)),lfp.data(ismember(lfp.timestamps,SWpeaks)),'r.')
    plot(lfp.timestamps(ismember(lfp.timestamps,SWpeaks)),lfp.data(ismember(lfp.timestamps,SWpeaks)),'g.')
    
        %OFF patches
%         y = get(gca,'Ylim');
%         patch([OFFints(:,1) OFFints(:,2) OFFints(:,2) OFFints(:,1)]',...
%             repmat([y(1) y(1) y(2) y(2)],length(OFFints(:,1)),1)',...
%             'r','FaceAlpha',0.1,'EdgeColor','none');
    
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
        %OFF patches
%         y = get(gca,'Ylim');
%         patch([OFFints(:,1) OFFints(:,2) OFFints(:,2) OFFints(:,1)]',...
%             repmat([y(1) y(1) y(2) y(2)],length(OFFints(:,1)),1)',...
%             'r','FaceAlpha',0.1,'EdgeColor','none');    
        %DOWN patches
        patch([DOWNints(:,1) DOWNints(:,2) DOWNints(:,2) DOWNints(:,1)]',...
            repmat([y(1) y(1) y(2) y(2)],length(DOWNints(:,1)),1)',...
            'g','FaceAlpha',0.2,'EdgeColor','none');
     
        
   % plot(deltaLFP.timestamps,normDELTA-peakthresh,'b','Linewidth',1)    
    plot(deltaLFP.timestamps(sampleIDX),normGAMMA(sampleIDX)+GAMMAwinthresh,'g','Linewidth',1) 
    plot(deltaLFP.timestamps(sampleIDX),normDELTA(sampleIDX)-DELTAwinthresh,'b','Linewidth',1) 
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

detectionparms.DELTApeakthresh = DELTApeakthresh;
detectionparms.DELTAwinthresh = DELTAwinthresh;
detectionparms.GAMMAdipthresh = GAMMAdipthresh;
detectionparms.GAMMAwinthresh = GAMMAwinthresh;


SlowWave.ints.UP = UPints;
SlowWave.ints.DOWN = DOWNints;
SlowWave.SWpeaks = SWpeaks;
SlowWave.SWpeakmag = SWpeakmag;
SlowWave.detectorinfo.detectorname = 'DetectSlowWaves';
SlowWave.detectorinfo.detectionparms = detectionparms;
SlowWave.detectiondate = today('datetime');

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