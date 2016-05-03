function [ spiketimes ] = JitterSpiketimes(spiketimes,jitterwin)
%[ spiketimes ] = JitterSpiketimes(spiketimes,jitterwin) jitters spiketimes
%of all cells within a designated window.
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] vectors of spiketimes for
%               each cell
%   jitterwin   time window within which to jitter spike times
%
%OUTPUT
%   spiketimes  jittered spiketimes
%
%
%DLevenstein 2016
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

