function [ stPR, stxcorr,t_corr ] = SpkTrigPopRate( spkObject,LFP,intervals,states,si_LFP )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% from Okun et al 2015
%
%
%TO DO:
%   -Add description and all that
%   -update code for ease of state vs not state.
%   -this is a memory suck..... fix it.  don't need to recalculate state
%   spikemat for each cell.... can split into a cell array earlier or
%   something
%   -smoothing seems to give different results from binning... why?
%   -LFP should be filtered... 0.1 - 200
%   -LFP should be downsampled to 200Hz?
%
%Last Updated: 7/23/15
%DLevenstein

%%

numcells = length(spkObject);


%Restrict spkObject,intervals to restrictioninterval
if exist('restrictinginterval','var')
    for i = states
        intervals{i} = intersect(intervals{i},restrictinginterval);
    end
end

if isobject(spkObject)
    for c = 1:numcells
        spiketimes{c} = Range(spkObject{c},'s');
    end
else
    spiketimes = spkObject;
end

%Downsample LFP... make this more flexible
downsample_factor = 10;
LFP = downsample(LFP,downsample_factor);
dt = si_LFP*downsample_factor;

display('Converting to Spikerate Matrix...')
T = [0 max(vertcat(spiketimes{:}))];
[spikemat,t] = SpktToSpkmat(spiketimes, T, dt);


%If want states
%Make Score vector
scorevec = zeros(size(t));
for s = states
    
    %For intervals in TSObjects
    statestarts = Start(intervals{s},'s'); 
    stateends = End(intervals{s},'s');
    
%     %For intervals not in TSObjects
%     statestarts = intervals(:,1);
%     stateends = intervals(:,2);
    
        for e = 1:length(statestarts)
            stateind = find(t>statestarts(e) & t<stateends(e));
            scorevec(stateind) = s;
        end
end


numstates = length(states);

tlag = 10;   %1s lag bounds for xcorr
tlag_dt = tlag/dt;
t_corr = [-tlag:dt:tlag];
stPR = zeros(numcells,numstates);
stxcorr = zeros(numcells,numstates,length(t_corr));
for c = 1:numcells
    display(['Calculating stLFP: Cell ',num2str(c),' of ',num2str(numcells)])
    for p = 1:numstates
        s = states(p); 
        
        %state specific spikes and LFP
        state_spikemat = spikemat(scorevec==s,:);   %spikemat for only times in state s
        stateLFP = LFP(scorevec==s);
        %Z-score the LFP in this state
        stateLFP = zscore(stateLFP);
        
        %Cross correlation - spikes with LFP, stPR is xcorr at 0 time lag
        stxcorr(c,p,:) = (1/sum(state_spikemat(:,c))).*xcorr(state_spikemat(:,c),stateLFP,tlag_dt);
        stPR(c,p) = stxcorr(c,p,t_corr==0);
    end
end

