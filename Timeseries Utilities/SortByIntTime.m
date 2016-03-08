function [intts,sortindex,intts_sorted,overlaps] = SortByIntTime(ts,ints,align)
%[intts,sortindex,intts_sorted] = SortByIntTime(timestamps,ints,align)
%sorts timestamps ts by their location relative to intervals ints.
%NOTE: as is, does not play nice with overlapping intervals
%
%INPUTS
%   timestamps
%   ints        [Nints x 2] array of interval onsets and offsets
%               (optional) can be TSObject
%   align       'onset', 'offset', or 'norm'
%
%OUTPUTS
%   intts       interval-relative timestamps
%   sortindex   index for sorting timestamps by interval-relative time
%   ints_sorted interval-relative timestamps, sorted
%   overlaps    cell array of other intts for those with overlaps
%
%
%TO DO:
%   -make play nicer with overlapping intervals
%
%Last Updated: 2/29/16
%DLevenstein
%%
intts = NaN(size(ts));
overlaps = cell(size(intts));

if isa(ints,'intervalSet')
    ints = [Start(ints,'s'), End(ints,'s')];
end

for ii = 1:length(ints(:,1))
    ininttimestamps = find(ts>=ints(ii,1) & ts<=ints(ii,2));
    newtimestamps = ts(ininttimestamps);
    
    if strcmp(align,'offset')
        newtimestamps = newtimestamps - ints(ii,2);
    elseif strcmp(align,'onset')
        newtimestamps = newtimestamps - ints(ii,1);
    elseif strcmp(align,'norm')
        intlength = ints(ii,2) - ints(ii,1);
        newtimestamps = (newtimestamps - ints(ii,1))/intlength;
    else
        display('var align must be "offset", "onset", or "norm"')
    end
    
    %If any of the timestamps are already in other intervals (overlapping)
    overlappingts = ininttimestamps(~isnan(intts(ininttimestamps)));
    if overlappingts
        if isempty(overlaps)
            display('Overlapping intervals...')
        end
        for tt = 1:length(overlappingts)
            overlaps{overlappingts(tt)} = [overlaps{overlappingts(tt)}, intts(overlappingts(tt))];
        end
    end
        
    
    intts(ininttimestamps) = newtimestamps;    
end

[intts_sorted,sortindex] = sort(intts);


end

