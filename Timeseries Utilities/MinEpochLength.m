function [ S, E, epoch_lengths ] = MinEpochLength( S, E, min_length, si)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%In case data starts in the middle of an epoch, drop first trigger
if E(1) < S(1)
    S = S(1:end-1);
    E = E(2:end);
else
end

% %In case data end in the middle of an epoch, drop last trigger
% if S(end) < E(end)
%     S = S(1:end-1);
%     E = E(2:end);
% else
% end

%How long are epochs?
epoch_length = E - S; 
epoch_length_t = epoch_length*si;

%Which epochs are longer than min_length
bigEpoch_indices = find(epoch_length_t > min_length);

S = S(bigEpoch_indices);
E = E(bigEpoch_indices);
epoch_lengths = E-S;

end

