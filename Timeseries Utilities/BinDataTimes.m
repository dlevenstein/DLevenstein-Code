function [ bincenters,binmeans,binstd,binnum ] = BinDataTimes( data,times,binedges,varargin )
%[ bincenters,binmeans,binstd,binnum ] = BinDataTimes( data,times,binedges )
%%
numbins = length(binedges)-1;
numdatacols = length(data(1,:));
%%
bincenters = binedges(1:end-1)+diff(binedges)/2;
binmeans = nan(numbins,numdatacols);
binstd = nan(numbins,numdatacols);
binnum = nan(numbins,numdatacols);
for bb = 1:numbins
    bintimes = times>=binedges(bb) & times<=binedges(bb+1);
    
    if ismember(varargin, 'sum')
        binmeans(bb,:) = sum(data(bintimes,:));
    else
    	binmeans(bb,:) = nanmean(data(bintimes,:));
    end
    
    binnum(bb,:) = numel(data(bintimes,:));
    binstd(bb,:) = nanstd(data(bintimes,:));
end


end

