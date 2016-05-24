function [ INT ] = IDXtoINT( IDX ,numstates)
%IDXtoINT(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    [t x 1] vector of state indices, where states are identified by
%           integers starting from 1.
%   numstates (optional)  number of interval types (for use
%
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times
%
%DLevenstein 2015-16
%%

if exist('numstates','var')
    states = 1:numstates;
else
    states = unique(IDX);
    numstates = length(states);
end

IDX = [0; IDX; 0];
for ss = 1:numstates
    statetimes = IDX==states(ss);
    INT{ss} = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
end



end

