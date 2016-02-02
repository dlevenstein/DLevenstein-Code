function [paircorrs,pairvalues] = PairMatToPairValues(corrmat,values)
%[paircorrs,pairvalues] = PairMatToPairValues(corrmat,values) 
%
%INPUT
%   corrmat     corrmat(i,j) is the weight of edge from element j to
%               element i
%   values      are a set of values for nodes
%
%Last Updated: 11/9/15
%DLevenstein
%%
numels = length(values);

%make sure values are vertical
if length(values(:,1))==1
    values = values';
end

ivaluemat = repmat(values,1,numels);
jvaluemat = repmat(values',numels,1);

paircorrs = corrmat(:);
pairvalues = [ivaluemat(:) jvaluemat(:)];

end

