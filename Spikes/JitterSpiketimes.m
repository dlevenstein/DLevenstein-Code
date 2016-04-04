function [ spiketimes ] = JitterSpiketimes(spiketimes,jitterwin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%%
numcells = length(spiketimes);
if isa(spiketimes,'tsdArray')
    for c = 1:numcells
        spiketimestemp{c} = Range(spiketimes{c},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
end


%Do this with cellfun....
for c = 1:numcells
    spiketimes{c} = spiketimes{c}+2*jitterwin*rand(size(spiketimes{c}))-jitterwin;
end



end

