function [NaNMat,t_align] = CellToNaNMat(epochcell,alignto,cut,sf)
%[NaNMat,t_align] = CellToNaNMat(epochcell,alignto,cut,sf) 
%   Detailed explanation goes here
%
%INPUT
%   epochcell   cell array of time series epochs
%   alignto     [Nepochs x 1] or [1 x 1] vector of alignment times (s).
%               negative number is align to that many from the end
%   cut         [Nepochs x 1] or [1 x 1] of amount of time to cut off the 
%               opposite end of the alignto time (i.e. to not include any
%               overhang window in PETH etc.)
%   sf          sampling frequency
%
%OUTPUT
%   NaNMat
%   t_align
%
%Last Updated: 10/6/15
%DLevenstein
%%

epochlengths = cellfun(@length,epochcell);
numepochs = length(epochcell);
maxlength = max(epochlengths);

if numepochs==0
    NaNMat = [];
    t_align = [];
    return
end

%For uniform alignment/cut times.
if length(alignto)==1
    alignto = repmat(alignto,1,numepochs);
end
alignto = alignto.*sf;
if length(cut)==1
    cut = repmat(cut,1,numepochs);
end
cut = cut.*sf;



NaNMat = NaN(2*maxlength,numepochs);
for ee = 1:numepochs
    %Align from end
    if alignto(ee)<0
        alignto(ee) = epochlengths(ee)+alignto(ee)-cut(ee);
        cellind = [(cut(ee)+1):epochlengths(ee)];
        matind = [1:(epochlengths(ee)-cut(ee))]+maxlength-alignto(ee);
    %Align from beginning    
    else
        cellind = [1:(epochlengths(ee)-cut(ee))];
        matind = [1:(epochlengths(ee)-cut(ee))]+maxlength-alignto(ee);
    end
    
    %Put it all in the NaNMat at the apropriate aligned location
    matind = int64(matind); %Bug fix...
    cellind = int64(cellind);
    NaNMat(matind,ee) = epochcell{ee}(cellind);
end

t_align = ([1:(2*maxlength)]-(maxlength))/sf;

nanrows = isnan(nanmean(NaNMat,2));
leftend = find(nanrows,1,'first')-1;
rightend = find(nanrows,1,'last')+1;
t_align([1:leftend,rightend:end]) = [];
NaNMat([1:leftend,rightend:end],:) = [];
end

    