function [thresh,cross,bihist] = BimodalThresh(bimodaldata,maxthresh)
%[thresh,cross,bihist] = BimodalThresh(bimodaldata) takes bimodal time
%series data, calculates the threshold between the modes (i.e. UP vs DOWN states),
%and returns the crossing times (i.e. UP/DOWN onset/offset times)
%
%INPUTS
%   bimodaldata     vector of bimodal data
%
%OUTPUS
%   thresh          threshold
%   cross
%       .upints     intervals of UP states [onsets offsets]
%       .downints   intervals of DOWN states [onsets offsets]
%   bihist
%       .bins       bin centers
%       .hist       counts
%
%DLevenstein Spring 2016
%%
if ~exist('maxthresh','var')
    maxthresh = inf;
end

%Remove data over the threshold... this is klugey
overmax = bimodaldata>=maxthresh;
bimodaldata(overmax) = nan;

numpeaks = 1;
numbins = 10; 
while numpeaks ~=2
    [bihist.hist,bihist.bins]= hist(bimodaldata,numbins);

    [PKS,LOCS] = findpeaks([0 bihist.hist 0],'NPeaks',2);
    LOCS = sort(LOCS)-1;
    numbins = numbins+1;
    numpeaks = length(LOCS);
    if numbins==25
    	cross.upints = []; cross.downints = []; thresh=nan; return
    end
end

betweenpeaks = bihist.bins(LOCS(1):LOCS(2));
[dip,diploc] = findpeaks(-bihist.hist(LOCS(1):LOCS(2)),'NPeaks',1,'SortStr','descend');

thresh = betweenpeaks(diploc);

%%
overind = bimodaldata>thresh;

crossup = find(diff(overind)==1);
crossdown = find(diff(overind)==-1);

%If there's only one crossing...
if isempty(crossup) || isempty(crossdown)
    cross.upints = []; cross.downints = []; return
end

%%
upforup = crossup;
upfordown = crossup;
downforup = crossdown;
downfordown = crossdown;

if crossup(1) < crossdown(1)
    upfordown(1) = [];
end
if crossdown(end) > crossup(end)
    downfordown(end) = [];
end
if crossdown(1) < crossup(1)
    downforup(1) = [];
end
if crossup(end) > crossdown(end)
    upforup(end) = [];
end
    


cross.upints = [upforup downforup];
cross.downints = [downfordown upfordown];


end

