function [popTE_in,popTE_out,...
    maxdelay_in,maxdelay_out] = PopulationTE(spikes,int,timerange)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   timerange (ms)
%
%See: Shimono and Beggs 2014, Ito et al 2011
%https://code.google.com/p/transfer-entropy-toolbox/
%
%%

if ~exist('timerange','var')
    timerange = [1 100];
end 
delays = [timerange(1):timerange(2)];


    
%Convert Spiketimes to a spike matrix with 1ms bins
dt = 0.001;
[spikemat,t] = SpktToSpkmat(spikes, [], dt);
popspikes = sum(spikemat,2);
spikemat = [popspikes,spikemat];

%Restrict to spikes in int (if applicable)
if exist('int','var')
    intidx = INTtoIDX(int,length(t),1/dt);
    intidx = logical(intidx);
    spikemat = spikemat(intidx,:);
end

%Convert Within-State Spikemat to asdf for TE function
asdf = SparseToASDF(spikemat',dt*1000);

%Calculate Transfer Entropy
i_order =3; %postsynaptic
j_order =3; %presynaptic
[peakTE,CI,allTE,maxdelay] = ASDFTE_parallel(asdf, delays,i_order, j_order);
popTE_in = peakTE(1,2:end);
popTE_out = peakTE(2:end,1);
maxdelay_in = maxdelay(1,2:end);
maxdelay_out = maxdelay(2:end,1);


end

