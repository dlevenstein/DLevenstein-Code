function [ data ] = IntToNaN(data,ints,sf,inv)
%[data] = IntToNaN(data,ints,sf) replaces data in ints with NaN
%
%INPUTS
%   inv     (optional) if 'inv', turn everything NOT in ints to NaN
%
%%
% if isa(ints,'intervalSet')
%     ints = [Start(ints,'s'), End(ints,'s')];
% end
% 
% ints = ints*sf;
% numints = length(ints(:,1));


% for ii = 1:numints
%     data(ints(ii,1):ints(ii,2),:) = NaN;
% end

intIDX = INTtoIDX(ints,length(data),sf);
intIDX = logical(intIDX);

if exist('inv','var')
    if strcmp(inv,'inv')
        data(~intIDX,:) = NaN;
    else
        display('inv must be "inv"');
    end
else
    data(intIDX,:) = NaN;
end

end

