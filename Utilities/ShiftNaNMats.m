function [ nanmat ] = ShiftNaNMats( nanmat )
%ShiftNaNMats Takes a left-aligned NaN-padded mat of intervals and
%right-aligns it
%
%
%TO DO:
%   -Make able to go the reverse way as well
%   -Make able to cope with vertical-time matrices
%
%Last Updated: 9/1/15
%DLevenstein


%How many NaNs in each interval row?
numnans = sum(isnan(nanmat),2);

numints = length(nanmat(:,1));
for ii = 1:numints
    nanmat(ii,:) = circshift(nanmat(ii,:),[0,numnans(ii)+1]);
end




end

