function [ ints ] = MergeSeparatedInts( ints,minseparation )
%UNTITLED3 Summary of this function goes here
%
%INPUT
%   ints            (n_ints x 2) matrix of interval start and end times
%   minseparation   merge ints separated by less than or equal to min
%%
intseparation = ints(2:end,1)-ints(1:end-1,2);
smallsep = find(intseparation<=minseparation);

for s = flipud(smallsep)'
    ints(s,2) = ints(s+1,2);
    ints(s+1,:) = [];
end

end

