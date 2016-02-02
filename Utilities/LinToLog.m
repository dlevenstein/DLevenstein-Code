function [ f_log,f_ind ] = LinToLog( f_lin,nfreqs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
f_log = logspace(log10(f_lin(1)),log10(f_lin(end)),nfreqs);

f_ind = zeros(size(f_log));
for ff = 1:nfreqs    
    [~,f_ind(ff)] = min(abs(f_log(ff)-f_lin));
end
f_log = f_lin(f_ind);

end

