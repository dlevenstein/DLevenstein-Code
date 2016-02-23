function [ ] = WatchRec()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
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

%% Load stuff for dev.
display('Loading From Database')
load('Database/BWData/BWRat19_032513/FilteredLFP.mat')
load('Database/BWData/BWRat19_032513/BWRat19_032513_LFPdata.mat')
load('Database/BWData/BWRat19_032513/BWRat19_032513_Intervals.mat')
load('Database/BWData/BWRat19_032513/BWRat19_032513_SSubtypes.mat')
load('Database/BWData/BWRat19_032513/BWRat19_032513_UPSpikeStatsE.mat')
load('Database/BWData/BWRat19_032513/BWRat19_032513_SpindleSpikeStats.mat')

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

vidtime = [800 2300];

%vidtime = [1400 2300];

%%
sf_LFP = 1250;
t_LFP = [1:length(data.LFP)]./sf_LFP;

LFP = data.LFP(t_LFP>=vidtime(1) & t_LFP<=vidtime(2));
t_LFP = t_LFP(t_LFP>=vidtime(1) & t_LFP<=vidtime(2));

display('Filtering 0.2-150 Hz...')
filtbounds = [0.2 150];
LFP = FiltNPhase(LFP, filtbounds, sf_LFP );

%% Spectrogram
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

REM = [Start(intervals{5},'s'), End(intervals{5},'s')];
SWS = [Start(intervals{3},'s'), End(intervals{3},'s')];
Wake = [Start(intervals{1},'s'), End(intervals{1},'s')];

%% UPs and Spindles
load('Database/BWData/BWRat19_032513/BWRat19_032513_UPSpikeStatsE.mat')
UPs = [isse.intstarts isse.intends];
load('Database/BWData/BWRat19_032513/BWRat19_032513_SpindleSpikeStats.mat')
SPs = [isse.intstarts isse.intends];
%% testfig


%%
LFP = zscore(LFP);


%% Figure
figure
for t1 = vidtime(1):framesize:vidtime(2)


%t1 = 0;
viewwin = [t1 t1+windowsize];
viewwin_LFP = [round(viewwin(1)*sf_LFP):round(viewwin(2)*sf_LFP)]+1;
viewwin_LFP = t_LFP>=viewwin(1) & t_LFP<=viewwin(2);

    subplot(18,1,1)

        plot(REM',ones(size(REM')),'g',...
            SWS',ones(size(SWS')),'b',...
            Wake',ones(size(Wake')),'r',...
            t1,1,'w.',...
            'LineWidth',10)
        ylim([0.99 1.01])
        xlim(vidtime)
        set(gca,'XTickLabel',[]);set(gca,'XTick',[]);
        set(gca,'YTickLabel',[])

    subplot(16,1,10:12)
        plot(REM', 2.4*ones(size(REM')),'g',...
            SWS', 2.4*ones(size(SWS')),'b',...
            Wake', 2.4*ones(size(Wake')),'r',...
            UPs(:,1),2.3*ones(size(UPs(:,1))),'+g',...
            UPs(:,2),2.3*ones(size(UPs(:,2))),'+r',...
            SPs(:,1),2*ones(size(SPs(:,1))),'+c',...
            SPs(:,2),2*ones(size(SPs(:,2))),'+m',...
            t_LFP(viewwin_LFP),LFP(viewwin_LFP),'k')

        xlim(viewwin)
        ylim([-3 3])

    
    subplot(16,1,2:9)
        imagesc(t_LFP(viewwin_LFP),log2(freqs),zspec(:,viewwin_LFP))
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        set(gca,'XTickLabel',[])
        xlim(viewwin)
        caxis([-1 4])
        axis xy

        
    subplot(4,1,4)
        plot(spindices(:,1),spindices(:,2),'.',...
            Ispindices(:,1),-Ispindices(:,2),'r.','MarkerSize',6)
        xlim(viewwin);ylim([-numIcells numcells+1]);
        xlabel('t (s)');
        ylabel('Neuron, Sorted by Rate')
        
        
        pause(pausetime)
end


end

