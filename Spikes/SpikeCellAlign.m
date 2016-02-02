function [ spikecell ] = SpikeCellAlign( spiketimes, eventtimes, windowlength )
% spikecell = SpikeCellAlign(spiketimes,eventtimes,windowlength) 
%
%
%Last Updated: 6/2/15
%DLevenstein


numevents = length(eventtimes);
numcells = length(spiketimes(1,:));


for e = 1:numevents
    for c = 1:numcells
        cellspikes = spiketimes{c};
        spikecell{e,c} = cellspikes(cellspikes>eventtimes(e) & cellspikes<eventtimes(e)+windowlength);
    end
end




end

