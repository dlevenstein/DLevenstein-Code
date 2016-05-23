function [ints,keepints,containingints] = RestrictInts(ints,restrictints)
%[ints,keepints,containingints] = RestrictInts(ints,restrictingints) 
%restricts intervals to those only appearing in restrictingints. 
%Can also return indices of restricting ints that contain ints
%
%INPUTS
%   ints            [Ni x 2] interval start and end times 
%                   or [Ni x 1] event times.
%   restrictints    [Nr x 2] interval start and end times
%
%OUTPUTS
%   ints            [No x 2] interval start and end times - only of
%                   intervals from ints that are within the intervals given
%                   by restrictints
%
%
%Last Updated: 9/21/15
%DLevenstein
%%
if isa(ints,'intervalSet')
    ints = [Start(ints,'s'), End(ints,'s')];
end
if isa(restrictints,'intervalSet')
    restrictints = [Start(restrictints,'s'), End(restrictints,'s')];
end

keepints = false(size(ints,1),1);
containingints = false(size(restrictints,1),1);

if isempty(ints) || isempty(restrictints)
    ints = [];
    return
end

if length(ints(1,:))==1
    ints = [ints ints];
end

%Loop over each restriction interval
numres = length(restrictints(:,1));
for rr = 1:numres
    resstart = restrictints(rr,1);
    resend = restrictints(rr,2);
    %Keep only intervals that start and end within a restriction interval
    keepints(ints(:,1)>=resstart & ints(:,2)<=resend) = true;
    if sum(ints(:,1)>=resstart & ints(:,2)<=resend)>=1
        containingints(rr) = true;
    end
end
ints = ints(keepints,:);
end

