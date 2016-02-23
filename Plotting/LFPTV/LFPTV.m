function [ output_args ] = LFPTV(datasetfolder,recname)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%

%Input
%   Spikes - (NCells x 1) Cell array of spiketimes
%   LFP
%   intervals - structure
%   events - structure
%   continuous variable(s)
%
%
%TO DO
%   -Add option to zoom (color, scale), go back, etc during playback, speed
%   up, slow down, pause
%   -turn into general purpose function with inputs as above
%   -add movement

%% DEV

recname = '20140526_277um';
% recname = 'Dino_061814_mPFC';
datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';



%% PrepLFPTV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([datasetfolder,recname,'/',recname,'_SSubtypes.mat'])
load([datasetfolder,recname,'/',recname,'_LFP.mat'])
load([datasetfolder,recname,'/',recname,'_StateIntervals.mat'])

%%
numcells = length(Se);
numIcells = length(Si);
for c = 1:numcells
    CellSpikes{c} = Range(Se{c},'s');
end
for c = 1:numIcells
    ICellSpikes{c} = Range(Si{c},'s');
end


%Sort by mean firing rate    
[~,sort_rate] = sort(cellfun(@numel,CellSpikes));
[~,Isort_rate] = sort(cellfun(@numel,ICellSpikes));
CellSpikes = CellSpikes(sort_rate);
ICellSpikes = ICellSpikes(Isort_rate);


%% Parms
windowsize = 10; % s
framesize = 0.01; %s
timedilation = 0.01; %not accurate... is dependent on CPU etc
pausetime = timedilation*framesize; %frames/s

%vidtime = [800 2300];
vidtime = [4315 4425];
%vidtime = [1400 2300];

%%
allLFP = LFP;

sf_LFP = 1250;
t_LFP = [1:length(LFP)]./sf_LFP;

LFP = LFP(t_LFP>=vidtime(1) & t_LFP<=vidtime(2));
t_LFP = t_LFP(t_LFP>=vidtime(1) & t_LFP<=vidtime(2));

display('Filtering 0.2-150 Hz...')
filtbounds = [0.2 150];
LFP = FiltNPhase(LFP, filtbounds, sf_LFP );

LFP = zscore(LFP);

%% Spectrogram

%Note: downsample to nernst for spectrogram
display('Calculating Wavelet Spectrogram...')
frange = [1 128];
nfreqs = 100;
ncyc = 10;
[ freqs, t_spec, spec ] = WaveSpec(LFP,frange,nfreqs,ncyc,1/sf_LFP,'log');

spec = abs(spec);
zspec = zscore(spec,[],2);

%% Get Spikeindices
display('Calculating Spikemat..')
[spikemat,t_spkmat,spindices] = SpktToSpkmat(CellSpikes, [0 t_LFP(end)], 0.1);

[~,~,Ispindices] = SpktToSpkmat(ICellSpikes, [0 t_LFP(end)], 0.1);

%% Intervals

REM = [Start(StateIntervals.REM,'s'), End(StateIntervals.REM,'s')];
SWS = [Start(StateIntervals.SWS,'s'), End(StateIntervals.SWS,'s')];
Wake = [Start(StateIntervals.Wake,'s'), End(StateIntervals.Wake,'s')];

%% UPs and Spindles

UPs = [Start(StateIntervals.UPstates,'s') End(StateIntervals.UPstates,'s')];
DOWNs = [Start(StateIntervals.DNstates,'s') End(StateIntervals.DNstates,'s')];
SPs = [Start(StateIntervals.Spindles,'s') End(StateIntervals.Spindles,'s')];

%% Global Spectrum

    %freqlist = linspace(0.5,55.5,111);
    window = 10;
    noverlap = 9;
    window = window*sf_LFP;
    noverlap = noverlap*sf_LFP;
    [FFTspec,FFTfreqs,t_FFT] = spectrogram(allLFP,window,noverlap,freqs,sf_LFP);
    FFTspec = abs(FFTspec);
    [zFFTspec,mu,sig] = zscore(log10(FFTspec)');





%% PlayLFPTV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
globalmap.spec 
globalmap.state



spikets
cellgroups

%% 
lfpmap.spec = zspec;
lfpmap.LFP = LFP;
lfpmap.freqs = freqs;
lfpmap.t_LFP = t_LFP;

spkidx.pE = spindices;
spkidx.pI = Ispindices;

stateints.REM = REM;
stateints.SWS = SWS;
stateints.Wake = Wake;

eventints.UP = UPs;
eventints.DOWN = DOWNs;
eventints.SP = SPs;

plotparms.statecolors = {'r','b','k'};
plotparms.eventcolors = {'g','r','m'};
plotparms.spkgroupcolors = {'k','r'};



tparms.vidtime = vidtime;
tparms.framedt = framesize;
tparms.framewidth = windowsize;

%%
statenames = fieldnames(stateints);
numstates = length(statenames);

eventnames = fieldnames(eventints);
numevents = length(eventnames);

cellgroupnames = fieldnames(spkidx); 
numcellgroups = length(cellgroupnames);

%%
%idea!  loop only over xlim....
%% Figure
figure
%for t1 = tparms.vidtime(1):tparms.framedt:tparms.vidtime(2)
t1 = tparms.vidtime(1);

viewwin = [t1 t1+tparms.framewidth];
viewwin_LFP = [round(viewwin(1)*sf_LFP):round(viewwin(2)*sf_LFP)]+1;
viewwin_LFP = lfpmap.t_LFP>=viewwin(1) & lfpmap.t_LFP<=viewwin(2);

%Global Recording Ref Map
    subplot(6,2,1)
        imagesc(t_FFT,log2(freqs),log10(FFTspec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)])
        %xlim(viewwin)
        %colorbar('east')
        ylim([log2(freqs(1)) log2(freqs(end))+0.2])
        ylabel({'LFP - FFT','f (Hz)'})
    
    hold on
    for ss = 1:numstates
        stateint = stateints.(statenames{ss});
        plot(stateint',log2(freqs(end))*ones(size(stateint')),plotparms.statecolors{ss},'LineWidth',8)
    end
        plot([t1 t1],get(gca,'ylim'),'w')
            
        %ylim([0.99 1.01])
        %xlim(vidtime)
        set(gca,'XTickLabel',[]);set(gca,'XTick',[]);
        set(gca,'YTickLabel',[])

%Wavelet Spectrogram    
    subplot(6,1,2:4)
        imagesc(lfpmap.t_LFP(viewwin_LFP),log2(lfpmap.freqs),lfpmap.spec(:,viewwin_LFP))
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        set(gca,'XTickLabel',[])
        xlim(viewwin)
        caxis([-1 4])
        axis xy

%LFP and Spikes     
    subplot(6,1,5:6)    
        hold on
    for ss = 1:numstates
        stateint = stateints.(statenames{ss});
        lineheight = 4;
        plot(stateint',lineheight*ones(size(stateint')),...
            plotparms.statecolors{ss},'LineWidth',3)
    end
    for ee = 1:numevents
        lineheight = 3.5-0.1*ee;
        eventint = eventints.(eventnames{ee});
        plot(eventint',lineheight*ones(size(eventint')),...
            plotparms.eventcolors{ee},'LineWidth',1)
    end
        plot(lfpmap.t_LFP(viewwin_LFP),lfpmap.LFP(viewwin_LFP),'k')

        xlim(viewwin)
    
    hold on
    groupoffset = -40;
    cellnumscale = 10;
    for cc = 1:numcellgroups
        groupspks = spkidx.(cellgroupnames{cc});
        plot(groupspks(:,1),(groupoffset-groupspks(:,2))./cellnumscale,...
            '.','Color',plotparms.spkgroupcolors{cc},'MarkerSize',8)
        groupoffset = groupoffset-max(groupspks(:,2));
    end
        
        xlim(viewwin);ylim([(groupoffset-1)/cellnumscale 4]);
        xlabel('t (s)');
        ylabel('Neuron, Sorted by Rate')
        
        
        pause(pausetime)
%end




end

