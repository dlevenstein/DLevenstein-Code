function [ stPR, stxcorr,t_corr ] = SpkTrigPopRate( spkObject,intervals,states,dt,restrictinginterval )
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
%
%Note: Removed input arg 'T' 7/14/15, will have to fix old scripts...
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
   

display('Converting to Spikerate Matrix...')



T = [0 max(vertcat(spiketimes{:}))];
[spikemat,t] = SpktToSpkmat(spiketimes, T, dt);

%Smoothed FR... matches binned very well in test dataset
% dt = 0.01;
% width =(12/sqrt(2))*0.001;% 0.001;
%[spikemat,t] = SpktToRate(spiketimes,width,T,dt);

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

tlag = 1;
tlag_dt = tlag/dt;
t_corr = [-tlag:dt:tlag];
stPR = zeros(numcells,numstates);
stxcorr = zeros(numcells,numstates,length(t_corr));
for c = 1:numcells
    display(['Calculating stPR: Cell ',num2str(c),' of ',num2str(numcells)])
    %Number of spikes for normalization, replaced with norm of spike mat
    %below
    numspikes = numel(spiketimes{c});   
    othercells = setdiff(1:numcells,c);
    %poprate = sum(spikemat_0mean(:,othercells),2);
    for p = 1:numstates
        s = states(p); 
        
        %state specific population rate
        state_spikemat = spikemat(scorevec==s,:);   %spikemat for only times in state s
        t_state = t(scorevec==s);   %t vector for state s
%         %center each cell around mean
 %        spikemat_0mean = state_spikemat-(ones(length(t_state),1)*mean(state_spikemat));
%         %sum all other cells
%         poprate = sum(spikemat_0mean(:,othercells),2);
        
        poprate = sum(state_spikemat(:,othercells),2);

        
       %try normalizing to pop rate mean
       %poprate = poprate./dt;
       %poprate = poprate-(ones(length(t_state),1)*mean(poprate));       
       %Zscore population rate
       poprate = zscore(poprate);
        
        stPR(c,p) = (1/sum(state_spikemat(:,c))).*(state_spikemat(:,c)'*poprate);
        stxcorr(c,p,:) = (1/sum(state_spikemat(:,c))).*xcorr(state_spikemat(:,c),poprate,tlag_dt);
    end
end

