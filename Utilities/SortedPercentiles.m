function [ datapercmean,datapercstd ] = SortedPercentiles(data,sortorder,numdistbins)
%SortedPercentiles(data,sort,numperc) averages data into equal number of
%bins as ordered by sortorder.
%
%
%
%Last Updated: 10/8/15
%DLevenstein
%%

numcells = length(data(1,:));

numperc = numcells/numdistbins;
percentilefloor = floor(linspace(1,numcells,numdistbins+1));
percentileceil = ceil(linspace(1,numcells,numdistbins+1));

numtbins = length(data(:,1));
datapercmean = zeros(numtbins,numdistbins);
datapercstd = zeros(numtbins,numdistbins);
data = data(:,sortorder);
for d = 1:numdistbins
    datapercmean(:,d) = nanmean(data(:,percentileceil(d):percentilefloor(d+1)),2);
    datapercstd(:,d) = nanstd(data(:,percentileceil(d):percentilefloor(d+1)),[],2);
end


end

