function [intts,sortindex,intts_sorted] = SortByIntTime(ts,ints,align)
%[intts,sortindex,intts_sorted] = SortByIntTime(timestamps,ints,align)
%sorts timestamps ts by their location relative to intervals ints.
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
%
%
%Last Updated: 11/24/15
%DLevenstein
%%
intts = NaN(size(ts));

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
    
    intts(ininttimestamps) = newtimestamps;    
end

[intts_sorted,sortindex] = sort(intts);


end

