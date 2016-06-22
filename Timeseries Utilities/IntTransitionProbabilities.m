function [ transprob ] = IntTransitionProbabilities(ints)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUTS
%   ints    {ninttypes} cell array of [nints x 2] set of start and end
%           times
%
%OUTPUTS
%   transprob   [ninttypes x ninttypes] matrix of probilities of
%               transitions. transprob(i,j) is the probability of
%               transition from state i to state j.
%
%DLevenstein Summer 2016
%%
numinttypes = length(ints);

transprob = zeros(numinttypes);
for ii = 1:numinttypes
    for jj = 1:numinttypes
        [~,~,ints1idx] = FindIntsNextToInts(ints{ii},ints{jj});
        numtrans = length(ints1idx);
        numints = length(ints{ii}(:,1));
        transprob(ii,jj) = numtrans./numints;
    end
end
        

end

