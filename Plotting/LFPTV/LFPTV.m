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
windowsize = 8; % s
timedilation = 0.01; %s/s    not accurate... is dependent on CPU etc

framesize = 0.01; %s
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





%% PlayLFPTV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




spikets
cellgroups

%% 
globalmap.spec = FFTspec;
globalmap.t = t_FFT;
globalmap.freqs = FFTfreqs;

lfpmap.spec = zspec;
lfpmap.LFP = LFP;
lfpmap.freqs = freqs;
lfpmap.t_LFP = t_LFP;
lfpmap.sf = sf_LFP;

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
tparms.pausetime = pausetime;


%%
PlayLFPTV(globalmap,lfpmap,spkidx,stateints,eventints,plotparms,tparms);

end

