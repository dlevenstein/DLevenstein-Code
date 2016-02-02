function [ISPCt,phaselag] = ISPCt(sig1phase,sig2phase,dt,win,sf_sig)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%   win     seconds
%%
if iscell(sig1phase)
    celllengths = cellfun(@length,sig1phase);
    sig1phase = vertcat(sig1phase{:});
end


siglen = length(sig1phase);
dt_sf = dt*sf_sig;
win_sf = win*sf_sig;

halfwin = win_sf/2;

tcenters = dt_sf:dt_sf:siglen;
numts = length(tcenters);
%%
ISPCt = zeros(numts,1);
phaselag = zeros(numts,1);
for t_i = 1:numts
    twin = (tcenters(t_i)-halfwin):(tcenters(t_i)+halfwin);
    twin(twin<=0) = [];
    twin(twin>siglen) = [];
    if strcmp(sig2phase,'diff')
        [ISPCt(t_i),phaselag(t_i)] = ISPC(sig1phase(twin),'diff');
    else
        [ISPCt(t_i),phaselag(t_i)] = ISPC(sig1phase(twin),sig2phase(twin));
    end
end




if exist('celllengths','var')
    ISPCt = mat2cell(ISPCt,nfreqs,celllengths);
end
%%
% figure
%     plot(IPSCt)
end

