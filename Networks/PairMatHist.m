function [corrbin,bincenters] = PairMatHist(paircorrs,pairvalues,numbins,binrange)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   paircorrs   [Npairs]
%   pairvalues  [Npairs x 2]
%   numbins     (or binedges)
%   binrange    (optional)
%
%Last Updated: 11/10/15
%DLevenstein
%%

if length(numbins)>1
    binedges = numbins;
    numbins = length(binedges)-1;
else
    if ~exist('binrange','var')
        binedges = linspace(min(pairvalues(~isinf(pairvalues))),max(pairvalues(~isinf(pairvalues))),numbins+1);
    else
        binedges = linspace(binrange(1),binrange(2),numbins+1);
    end
end

corrbin.mean = zeros(numbins,numbins);
corrbin.std = zeros(numbins,numbins);
corrbin.num = zeros(numbins,numbins);

%%
ival = pairvalues(:,1);
jval = pairvalues(:,2);

bincenters = diff(binedges)/2+binedges(1:end-1);
binedges(1) = -Inf;
binedges(end) = Inf;
for f_i = 1:numbins
    for f_j = 1:numbins
        corrs = paircorrs(ival>=binedges(f_i) & ival<=binedges(f_i+1) & jval>=binedges(f_j) & jval<=binedges(f_j+1));
        corrbin.mean(f_i,f_j) = nanmean(corrs);
        corrbin.std(f_i,f_j) = nanstd(corrs);
        corrbin.num(f_i,f_j) = numel(corrs(~isnan(corrs)));
    end
end

corrbin.binedges = binedges;
corrbin.bincenters = bincenters;
end

