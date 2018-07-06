function [ bincenters,binmeans,binstd,binnum,binneddata ] = BinDataTimes( data,times,binedges,varargin )
%[ bincenters,binmeans,binstd,binnum,binneddata ] = BinDataTimes( data,times,binedges )
%%
if isrow(data); data = data'; end

if isrow(binedges) || iscolumn(binedges)
    binedges_temp(:,1) = binedges(1:end-1);
    binedges_temp(:,2) = binedges(2:end);
    binedges = binedges_temp;
end
bincenters = mean(binedges,2);

if isempty(data)
    binmeans = nan(size(bincenters));
    binstd=nan(size(bincenters));
    binnum=zeros(size(bincenters)); binneddata=[];
    return
end

numbins = length(binedges(:,1));
numdatacols = length(data(1,:));
%%

binmeans = nan(numbins,numdatacols);
binstd = nan(numbins,numdatacols);
binnum = nan(numbins,numdatacols);
for bb = 1:numbins
    bintimes = times>=binedges(bb,1) & times<=binedges(bb,2);
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

