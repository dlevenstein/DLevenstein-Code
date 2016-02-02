function [ INT ] = IDXtoINT( IDX )
%IDXtoINT(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    [t x 1] vector of state indices, where states are identified by
%           integers starting from 1.
%
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times



states = unique(IDX);
numstates = length(states);

IDX = [0; IDX; 0];
for ss = 1:numstates
    statetimes = IDX==states(ss);
    INT{ss} = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
end



end

