function [ bincenters,binmeans,binstd,binnum,binneddata ] = BinDataTimes( data,times,binedges,varargin )
%[ bincenters,binmeans,binstd,binnum,binneddata ] = BinDataTimes( data,times,binedges )
%%
if isrow(data); data = data'; end

numbins = length(binedges)-1;
numdatacols = length(data(1,:));
%%
bincenters = binedges(1:end-1)+diff(binedges)/2;
binmeans = nan(numbins,numdatacols);
binstd = nan(numbins,numdatacols);
binnum = nan(numbins,numdatacols);
for bb = 1:numbins
    bintimes = times>=binedges(bb) & times<=binedges(bb+1);
    binneddata{bb} = data(bintimes,:);
    
    if ismember(varargin, 'sum')
        binmeans(bb,:) = sum(binneddata{bb});
    else
    	binmeans(bb,:) = nanmean(binneddata{bb});
    end
    
    binnum(bb,:) = numel(binneddata{bb});
    binstd(bb,:) = nanstd(binneddata{bb});
end


end

