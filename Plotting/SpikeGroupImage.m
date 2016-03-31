function [ spikegroupimage ] = SpikeGroupImage(sitevalues,spikegroups)
%SiteGroupImage plots the values for each recording site in sitevalues as
%arranged in spikegroups.
%
%%
numgroups = length(spikegroups);
longestgroup = max(cellfun(@length,spikegroups));
spikegroupimage = nan(longestgroup,numgroups);
sitematidx = nan(longestgroup,numgroups);

for ss = 1:numgroups
    grouplen = length(spikegroups{ss});
    sitematidx(1:grouplen,ss) = spikegroups{ss}';
end
spikegroupimage = sitevalues(sitematidx+1);

imagesc(spikegroupimage)

end

