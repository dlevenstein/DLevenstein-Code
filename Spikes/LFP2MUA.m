function [ sp_indices, spikes, MUAPower, MUA, srate] = LFP2MUA(data,sampfreq)
%[sp_indices,spikes,MUAPower,MUA,srate] = LFP2MUA(LFP, sampfreq) 
%
%
%TO DO:
%   -Verify/fix gaussian normalization
%
%Last Updated: 3/31/15
%DLevenstein

SHOWFIGS = false;

%High-pass filter
spikepass = [400 7000];
[hpData,MUAPower] = FiltNPhase(data,spikepass,sampfreq);
hpData = -hpData;   %Flip high-passed LFP for spike identification
MUA = hpData;

%Z-score the highpassed spiketrace
hpData = (hpData-mean(hpData))/std(hpData);

%Threshold to Isolate Spike Times
th = 3.5;   %Threshold, (stds)
thresh_filter = hpData > th;
spikes = thresh_filter(2:end) > thresh_filter(1:end-1);
spikes = [0 ; spikes];          %vector of 0s 1s
sp_indices = find(spikes == 1); %indices of spikes

si = 1/sampfreq;
t_si = [0:length(data)-1];
t = t_si*si;

%Gaussian window
width = 0.01; %s
wind_x = [-2.5*width:si:2.5*width];
window = Gauss(wind_x,0,width);
window = (window/sum(window))/si;

srate = FConv(window,spikes');


%%Plot filtered and unfiltered spikes
if SHOWFIGS == true
    figure; numplots = 3;
        subplot(numplots,1,1)
            plot(t,data)
            ylabel('mV')
            xlim([56.4 56.6]);
        subplot(numplots,1,2)
            plot(t,15*spikes,'r', t,hpData,'k');
            xlim([56.4 56.6]);
end

end

