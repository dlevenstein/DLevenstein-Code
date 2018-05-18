function [ borders ] = ConditionBorders( metrics,bounds )
%Finds overlap in condition bounds for multiple metrics
%
%INPUTS
%   metrics     [Nconditions x [Ndim]]
%   bounds      [Nconditions x 2]
%
%
%DLevenstein
%%

numconditions = length(bounds(:,1));
for cc = 1:numconditions
    inbounds(cc,:,:) = metrics(cc,:,:)>=bounds(cc,1) & metrics(cc,:,:)<=bounds(cc,2);
end
   
fitconditions = squeeze(all(inbounds,1));

borders = bwboundaries(fitconditions);

% figure
% imagesc(fitconditions)
% hold on
% for k = 1:length(borders)
%    boundary = borders{k};
%    plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 1)
% end
% axis xy

end

