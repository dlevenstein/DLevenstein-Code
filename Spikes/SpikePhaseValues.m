function [spikephases,spikepowers] = SpikePhaseValues(spiketimes,LFP,frange,sf)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%
%
%INPUTS
%   spiketimes      {numspikes} cell array - each cell is spiketimes of a
%                   neuron. (Optional: can be tsdArray)
%
%
%Last Updated: 10/9/15
%DLevenstein
%%
numcells = length(spiketimes);

if isa(spiketimes,'tsdArray')
    for c = 1:numcells
        spiketimestemp{c} = Range(spiketimes{c},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
end
%%
[~,power,phase] = FiltNPhase(LFP,frange,sf);
power = zscore(power);
%%
%Convert spiketimes to samples
spikesamps = cellfun(@(X) round(X*sf),spiketimes,'UniformOutput',false);

%Find power and phase values at each spike
spikephases = cellfun(@(X) phase(X),spikesamps,'UniformOutput',false);
spikepowers = cellfun(@(X) power(X),spikesamps,'UniformOutput',false);



end

