function [stSpec,stSpec_raw,stSpec_std,stSpec_mean,freqs] = SpkTrigSpec2( spkObject,LFP,si_LFP,intervals,states,restrictinginterval)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Spike Triggered Spectrogram with Spike Shuffling
%
%
%INPUT
%   spkObject   spike times, cell array of spike times or TSObject
%   LFP         LFP signal
%   si_LFP      Sampling interval of the LFP - recommended: 0.001 (1kHz)
%
% (optional)
%   intervals   [n x 2] matrix of state start and end times. column 1 is
%               start times, column 2 is end times.
%               Alternatively, intervals can be TSObject
%   states      if intervals is TSObject, for multiple states. vector of
%               state ID numbers to calculate stSpec for
%   restrictinginterval     TSObject to limit to subset of time points
%
%
%OUTPUT
%   stSpec      [Ncells x Nfreqs x Nstates] spike triggered spectrum for 
%               each cell in each state, in units of standard deviation
%               relative to spike shuffled data
%   stSpec_raw  [Ncells x Nfreqs x Nstates] raw spike triggered spectrum
%               for each cell in each state
%   stSpec_std  [Ncells x Nfreqs x Nstates] std of spike-shuffled data
%               for each cell 
%   stSpec_mean [Nstates x Nfreqs] mean spectrum of spike-shuffled data for
%               each cell in each state
%   freqs       vector of frequencies
%
%
%Dependencies
%   SpktToSpkmat
%   WaveSpec
%   WaveFilt
%   MorletWavelet
%   FConv
%
%TO DO:
%   -Add description and all that
%   -might want to clean up state vs non-state...
%   -option to select freq range etc
%
%Last Updated: 9/10/15
%DLevenstein
%%
% Dev... BWData

% load('Database/BWData/BWRat19_032513/BWRat19_032513_SSubtypes.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_GoodSleepInterval.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_Intervals.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_LFPdata.mat')
% 
% 
% spkObject = Se;
% states = [1,3,5];
% 
% restrictinginterval = GoodSleepInterval;
% spkObject = Se;
% LFP = data.LFP;
% si_LFP = 1/1250;


%%
% % Dev... WMData
% load('Database/struct/15111004.mat')
% spkObject = {recdata.spikes};
% LFP = recdata.LFP;
% si_LFP = recdata.si;
% 
% intervals = [recdata.wh_on' recdata.wh_off'];
% intervals = [recdata.wh_off(1:end-1)' recdata.wh_on(2:end)'];

%%
display('Spike Triggered LFP Spectrum');
numcells = length(spkObject);
t_LFP = ([1:length(LFP)])*si_LFP;

%Restrict intervals to restrictioninterval
if exist('restrictinginterval','var')
    goodint = [Start(restrictinginterval,'s') End(restrictinginterval,'s')];
    for i = states
        intervals{i} = intersect(intervals{i},restrictinginterval);
    end
else
    goodint = [0 length(LFP)*si_LFP];
end

%Extract Spikes from TSObject to cell array 
%(N/A if spikes are already a cell array)
if isobject(spkObject)
    for c = 1:numcells
        spiketimes{c} = Range(spkObject{c},'s');
    end
else
    spiketimes = spkObject;
end

%If no designated states and/or intervals
if ~exist('states','var')
    states = 1;
end

if ~exist('intervals','var')
    intervals = [0 length(LFP)*si_LFP];
end




%%
LFP = LFP(t_LFP>=goodint(1) & t_LFP<=goodint(2));
t_LFP = t_LFP(t_LFP>=goodint(1) & t_LFP<=goodint(2));

% %Downsample LFP... make this more flexible. want LFP and spikemat to have
% %same dt?  or maybe just do this via spiketimes for now just accept
% si_LFP as dt... downsample for faster analysis and less RAM
% consumption... as long as sf_LFP is >=1kHz should be fine
downsample_factor = 4;
LFP = downsample(LFP,downsample_factor);
t_LFP = downsample(t_LFP,downsample_factor);
dt = si_LFP*downsample_factor;


%Make Spike Matrix
display('Converting to Spikerate Matrix...')
T = [goodint(1) goodint];
[spikemat,t] = SpktToSpkmat(spiketimes, T, dt);



%% Wavelet Transform LFP
display('Calculating Wavelet Spectrogram...')


%Check alignment t_LFP == t.  May want them to be same orientation?
%This alignment thing is a mess.....
% assert(t(1)==t_LFP(1) & t(end)==t_LFP(end) & length(t)==length(t_LFP),...
%     'LFP and Spiketimes misaligned... squish that bug...')

frange = [0.5 128];
nfreqs = 100;
ncyc = 10;
[freqs,~,LFPspec] = WaveSpec(LFP,frange,nfreqs,ncyc,dt,'log');
LFPspec = log10(abs(LFPspec));



%%
%If want states
%Make Score vector
scorevec = zeros(size(t));
for s = states   
    %For intervals in TSObjects
    if iscell(intervals)
        statestarts = Start(intervals{s},'s'); 
        stateends = End(intervals{s},'s');
    else
    %For intervals not in TSObjects
        statestarts = intervals(:,1);
        stateends = intervals(:,2);  
    end
        for e = 1:length(statestarts)
            stateind = find(t>statestarts(e) & t<stateends(e));
            scorevec(stateind) = s;
        end
end


numstates = length(states);



%Initiate Outputs
nshuff = 1000;
stSpec = zeros(numcells,nfreqs,numstates);
stSpec_raw = zeros(numcells,nfreqs,numstates);
stSpec_std = zeros(numcells,nfreqs,numstates);
stSpec_mean = zeros(numcells,nfreqs,numstates);
%In each state, get the spectrum at spike times
for p = 1:numstates
    s = states(p);
    state_spikemat = spikemat(scorevec==s,:);   %spikemat for only times in state s
    state_spec = LFPspec(:,scorevec==s);    %dimensions are f'd
    
    for c = 1:numcells
        %Spec at each spike time in state - mean and std
        cellspikespec = state_spec(:,find(state_spikemat(:,c))); 
        stSpec(c,:,p) = mean(cellspikespec,2);
        
        %Make circularly shuffled spikes
        nshuff = 1000;
        shuff = randi(length(state_spikemat(:,c)),[length(state_spikemat(:,c)),1]);
        stSpec_shuff = zeros(nfreqs,nshuff);
        for ns = 1:nshuff
            shuffspikes =  circshift(state_spikemat(:,c),shuff(ns));
        	cellspikespec = state_spec(:,find(shuffspikes)); 
            stSpec_shuff(:,ns) = mean(cellspikespec,2);
        end
        %Mean and std of shuffled spike triggered Spectrum
        stSpec_mean(c,:,p) = mean(stSpec_shuff,2);
        stSpec_std(c,:,p) = std(stSpec_shuff,[],2);
        stSpec_raw(c,:,p) = stSpec(c,:,p);
        stSpec(c,:,p) = (stSpec(c,:,p)-stSpec_mean(c,:,p))./stSpec_std(c,:,p);

    end

end







end

