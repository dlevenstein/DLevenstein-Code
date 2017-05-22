function [ symmetricweightmat ] = MakeSymmetric( weightmat )
%[ symmetricweightmat ] = MakeSymmetric( weightmat ) takes a directed
%weight matrix and makes it symmetric/undirected by adding the two corners
%together.
%
%
%DLevenstein 2017
%%

symmetricweightmat = weightmat+triu(weightmat,1)'+tril(weightmat,-1)';

%%

% figure
% subplot(2,2,1)
% imagesc(weightmat)
% subplot(2,2,2)
% imagesc(symmetricweightmat)
end

