function [intlagmean,lags] = XCorrInInt(signal,int,maxlag,sf)
%[intxcorr,lags] = XCorrInInt(signal,int,maxlag,sf) takes the
%cross-correlation of a signal with reference to only the time points in
%interval set int.
%
%
%DLevenstein Fall 2016
%%
% signal = spPOWER;
% sf = sf_best;
% int = IDXtoINT(signal>=1.5 & signal<=2);
% int = int{1}./sf;
% maxlag = 5;

%%
maxlag_si = round(maxlag.*sf);
lags = [-maxlag_si:maxlag_si]';

if isequal(size(signal),size(int))
    %Convert Logicals to Indices
    idx = find(int);
else
    %Convert Intervals to Indices
    idx = INTtoIDX({int},length(signal),sf);
    idx = find(idx);
end
%%
intlagmean = zeros(size(lags));
intxcorr = zeros(size(lags));

nanholder = length(signal)+1;
signal(nanholder) = nan;
for ii = 1:length(idx)
    windices = idx(ii) + lags;
    windices(windices<1 | windices>nanholder) = nanholder;
    intlagmean = nansum([intlagmean,signal(windices)],2);
    
    %An attempt here at doing a (more rigorous?) xcorr - this needs to be
    %verified etc etc.
    xcorrtemp = signal(windices).* signal(idx(ii));
    intxcorr = nansum([intxcorr,xcorrtemp],2);
end
intlagmean = intlagmean./length(idx);
lags = lags./sf;
%%
% figure
% plot(lags,intlagmean,'k')
end

