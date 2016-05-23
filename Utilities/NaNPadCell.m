function [ cellarray_padded ] = NaNPadCell(cellarray,side)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%

%%
epochlens = cellfun(@(X) length(X(:,1)),cellarray,'UniformOutput',false);
maxlen = max([epochlens{:}]);

lendiff = cellfun(@(X) maxlen-X,epochlens,'UniformOutput',false);

numcells = length(cellarray);
for cc = 1:numcells
    numcols = length(cellarray{cc}(1,:));
    if strcmp(side,'end')
        cellarray_padded{cc} = vertcat(cellarray{cc},NaN(lendiff{cc},numcols));
    elseif strcmp(side,'start')
        cellarray_padded{cc} = vertcat(NaN(lendiff{cc},numcols),cellarray{cc});
    else
        display('Which side do you want to put the NaNs on?')
        continue
    end
end

