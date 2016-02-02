function [ epochintmean ] = EpochIntMean(epochcell,int,sf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%%
int = sf.*int;
if min(int(:))<=0;
    display('Intervals less than 1, replacing with 1')
    int(find(int<1)) = 1;
end
%%
starts = num2cell(int(:,1));
ends = num2cell(int(:,2));

epochintmean = cellfun(@(A,B,C) mean(A(B:C)),epochcell,starts,ends,...
    'UniformOutput',false);
epochintmean =cell2mat(epochintmean);

end

