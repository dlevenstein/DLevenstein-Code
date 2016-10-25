function [allintcorr,intcorr] = CompareIntervalSets(intset1,intset2,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   intset1     {nstates} cell array of state intervals - start and end times
%   intset2     cell array of equal number of interval states
%   OPTIONAL
%       'nulltimes'     0 or 1 to not count or count unscored times
%                       (default is 1)
%       intsets can alternatively be intervalSet objects from TSObjects
%       toolbox
%
%
%
%OUTPUT
%   allintcorr  total correlation between intervals
%   intcorr     [nstates] vector of correlations for each state
%   mismatch    mismatching times
%
%
%IMPORTANT NOTE: ISSUE WITH OVERLAPPING INTERVALS
%
%Last Updated: 9/23/15
%DLevenstein
%%
if isa(intset1{1},'intervalSet')
    for ii = 1:length(intset1)
        intset1{ii} = [Start(intset1{ii},'s'), End(intset1{ii},'s')];
    end
end

if isa(intset2{1},'intervalSet')
    for ii = 1:length(intset2)
        intset2{ii} = [Start(intset2{ii},'s'), End(intset2{ii},'s')];
    end
end

if find(strcmp(varargin,'nulltimes'))
    nullarg = find(strcmp(varargin,'nulltimes'));
	nulltimes = varargin{nullarg+1};
else
    nulltimes = 1;
end


%%
allinttimes = vertcat(intset1{:},intset2{:});
idxlength = max(allinttimes(:));

idx1 = INTtoIDX(intset1,idxlength);
idx2 = INTtoIDX(intset2,idxlength);

if ~nulltimes
    zerotimes = (idx2==0 | idx1==0);
    idx1(zerotimes) = [];
    idx2(zerotimes) = [];
end

sameidx = idx1==idx2;
numsame = sum(sameidx);
allintcorr = numsame/length(sameidx);

intcorr = zeros(size(intset1));
for ss = 1:length(intset1)
    intidx1 = idx1==ss;
    intidx2 = idx2==ss;
    intidxeither = (intidx1==1 | intidx2==1);
    sameidx = (intidx1(intidxeither))==(intidx2(intidxeither));
    numsame = sum(sameidx);
    intcorr(ss) = numsame/length(sameidx);
end

%%
% figure
%     subplot(2,1,1)
%         plot(idx1)
%     subplot(2,1,2)
%         plot(idx2)




end

