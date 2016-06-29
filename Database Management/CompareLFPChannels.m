function [noisescore,noiseresults] = CompareLFPChannels(pLFP)
%[numbadseconds,noisetimes] = CompareLFPChannels(pLFP) detects different 
%types noise on multiple LFP channels.
%
%INPUT
%   pLFP [t x numchannels] vector of LFP timeseries
%
%OUTPUT
%   noisescore [numchannels] vector giving a rating for each pLFP column.
%               note: noisescore is internally consistent, cannot compare
%               across recordings. larger score = more noisy
%   noiseresults    structure giving details... need to add separate output
%                   for each kind of noise...
%       .numbadseconds
%       .noisetimes
%
%DLevenstein 2016
%%
absLFPthresh = 10; %Threshold for LFP magnitude. Units: standard Deviations
dtLFPthresh = 5; %Threshold for dLFP/dt magnitude. Units: standard Deviations/timebin

numchannels = length(pLFP(1,:));

pLFP = zscore(pLFP);

absLFP = abs(pLFP);
dtLFP = abs(diff(pLFP));

largeLFPnoise = absLFP>absLFPthresh;
LFPdtnoise = dtLFP>dtLFPthresh;% | dtLFP==0;   should only be if 0 for multiple time steps
LFPdtnoise = [LFPdtnoise; zeros(1,numchannels)];

noiseresults.noisetimes = largeLFPnoise | LFPdtnoise;
noiseresults.numbadseconds = sum(noiseresults.noisetimes,1);
noisescore = noiseresults.numbadseconds;

%%
% figure
%     subplot(4,1,1)
%         plot(pLFP)
%     subplot(4,1,2)
%         plot(dtLFP)
%     subplot(2,2,3)
%     hist(absLFP)
%     subplot(2,2,4)
%     hist(dtLFP)
    

end

