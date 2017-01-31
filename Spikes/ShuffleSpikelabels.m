function [ spiketimes ] = ShuffleSpikelabels(spiketimes)
%[ spiketimes ] = ShuffleSpikelabels(spiketimes) jitters spiketimes
%of all cells within a designated window.
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] vectors of spiketimes for
%               each cell
%
%OUTPUT
%   spiketimes  jittered spiketimes
%
%
%DLevenstein 2017
%%
numcells = length(spiketimes);
if isa(spiketimes,'tsdArray')
    for cc = 1:numcells
        spiketimestemp{cc} = Range(spiketimes{cc},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
end

%Count how many spikes each cell should get at the end
numcellspks = cellfun(@length,spiketimes);

%Put all the spike times together into a single vector, then shuffle them
allspks = cat(1,spiketimes{:});
allspks = allspks(randperm(length(allspks)));

%Put the shuffled times back into a cell array with same number of spikes
%per cell
spiketimes = mat2cell(allspks,numcellspks)';
spiketimes = cellfun(@sort,spiketimes,'UniformOutput',false);


end

