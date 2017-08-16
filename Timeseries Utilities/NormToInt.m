function [normdata,intmean,intstd] = NormToInt(data,normtype,int,sf,varargin)
%NormToInt(data,normtype, int,sf) normalizes the data to a subset of time 
%(int). Has options for multiple normalization types.
%
%INPUTS
%   data    [Nt x Ndim] 
%   normtype 'mean', 'Z', 'max', 'percentile', 'modZ' for modified Z-score
%   int     [Nints x 2] reference interval onsets and offset times   
%           to normalize the data with respect to (optional)
%   sf      (optional) sampling frequency of the data. default 1
%
%Note: modified Z score is median-based and robust to outliers
%see https://ibm.co/2qi4Vy5
%
%DLevenstein Summer 2016
%TO DO: Improve input parsing, compadible with buzcode structures
%%
SHOWFIG = false;


%%
if ~exist('int','var') || isempty(int)
    int = [1 size(data,1)];
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


numints = length(int(:,1));
if ~exist('sf','var')|| isempty(sf)
    sf = 1;
end


int = round(int*sf);
int(int==0)=1; %Turn time=0 to the first indexed datapoint

int_data = [];
for ii = 1:numints
    int_data = [int_data; data(int(ii,1):int(ii,2),:)];
end
%Instead turn non-int data to nan and take mean, can then use movmean etc
%for moving stuff...

intmean = nanmean(int_data);
intstd = nanstd(int_data);

intmedian = median(int_data);
intMAD = mad(int_data,1);

switch normtype
    case 'Z'
        normdata = (data-repmat(intmean,length(data(:,1)),1))./repmat(intstd,length(data(:,1)),1);
    case 'mean'
        normdata = bsxfun(@(X,Y) X./Y,data,intmean);
    case 'max'
        colmax = max(int_data,[],1);
        normdata = bsxfun(@(X,Y) X./Y,data,colmax);
    case 'percentile'
        sortdata = unique(sort(int_data(~isnan(int_data))));
        percentiles = linspace(0,1,length(sortdata));
        normdata = interp1(sortdata,percentiles,data,'nearest');
    case 'modZ'
        normdata = 0.6745*(data-repmat(intmedian,length(data(:,1)),1))./repmat(intMAD,length(data(:,1)),1);
        intmean = intmedian;
        intstd = intMAD;
    otherwise
        display('incorrect normtype')
end

%%
if SHOWFIG
    figure
        subplot(2,1,1)
            plot(data,normdata,'.')
        subplot(2,1,2)
            hist(int_data)
end


end

