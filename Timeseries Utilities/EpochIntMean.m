function [ epochintmean ] = EpochIntMean(epochcell,int,sf)
%[epochintmean] = EpochIntMean(epochcell,int,sf) takes a cell array of
%isolated epochs and finds the mean value in some interval in those epochs.
%
%%
numepochs = length(epochcell);
if length(int(:,1))==1
    int = repmat(int,numepochs,1);
elseif isempty(int)
    display('No Interval...')
    epochintmean = [];
    return
end

int = round(sf.*int);
if min(int(:))<=0;
    display('Intervals less than 1, replacing with 1')
    int(find(int<1)) = 1;
end

if isempty(epochcell)
    display('No Epochs...')
    epochintmean = [];
    return
end


%%
starts = num2cell(int(:,1));
ends = num2cell(int(:,2));

epochintmean = cellfun(@(A,B,C) mean(A(B:C,:)),epochcell,starts,ends,...
    'UniformOutput',false);
epochintmean =cell2mat(epochintmean);

end

