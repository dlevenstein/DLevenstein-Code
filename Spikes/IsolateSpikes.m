function [ sp_indices, spikes, t_wave, spikewaves, pt] = IsolateSpikes( data, sampfreq,varargin )
%IsoalteSpikes(data, sf) takes a cell attached and returns
%spike times as indices and as a 0-vector with 1s at spiketimes.
%
%Optional Arguement: Save As, thresh.   IN THAT ORDER!
%
%TO DO: 
%   fix issue with noisy input (see Cell 26,recording2)
%   -allow for arbitrarty threshold without saving
%   -include refractory period in figure
%
%Last Updated: 4/10/15
%DLevenstein

si = 1/sampfreq;
t_si = [0:length(data)-1];
t = t_si*si;


twindow = [0.001 0.002]; %window before and after (ms)
twindow_si = twindow/si;


if nargin == 4
    th = varargin{2};
else
    th = 9;   %Threshold, (stds)
end    
    
%High-pass filter
%spikepass = [200 9000]; %original
spikepass = [400 5000];
%[800-5000] fujisawa
filterorder = 4; %cycles lowest freq
[hpData] = FiltNPhase(data,spikepass,sampfreq,filterorder);
%Z-score the highpassed spiketrace
hpData = (hpData-mean(hpData))/std(hpData);

%Threshold to Isolate Spike Times
thresh_filter = hpData > th;
spikes = thresh_filter(2:end) > thresh_filter(1:end-1);
spikes = [0 ; spikes];          %vector of 0s 1s
sp_indices = find(spikes == 1); %indices of threshold crossings
%No Spikes with windows that hang over edge
numshort = sum(sp_indices-twindow_si(1) < 1);
numlong = sum(sp_indices+twindow_si(2) > length(hpData));
sp_indices = sp_indices((1+numshort):(end-numlong));
%No spikes within the window refractory period
ISIs = (sp_indices(2:end)-sp_indices(1:end-1))*si;
smallISIs = find(ISIs<twindow(2))+1;
sp_indices(smallISIs) = [];

numspikes = length(sp_indices);

if numspikes == 0
    display('No events above threshold')
    sp_indices = [];
    spikes = [];
    t_wave = [];
    spikewaves = [];
    pt = [];
    return
end


%Find Max within a window
spikewindow = ones(twindow_si(2),numspikes);
spikewindow(1,:) = sp_indices;
spikewindow = cumsum(spikewindow,1);
spikewaves = hpData(spikewindow);
[maxvalue, maxindex] = max(spikewaves,[],1);
[minvalue, minindex] = min(spikewaves,[],1);

%Remove spikes that are too large - most likely noise (make this
%optional...)
maxlimit = 100; %STD
bigspikes = find(maxvalue>maxlimit);
maxvalue(bigspikes) = [];
maxindex(bigspikes) = [];
minvalue(bigspikes) = [];
minindex(bigspikes) = [];


%Peak To Trough
pt = (minindex-maxindex)*si;
pt_sd = std(pt);
pt = mean(pt);

%Adjust Spike indices to peak
sp_indices = sp_indices' + maxindex;
spikewindow = ones(sum(twindow_si)+1,numspikes);
spikewindow(1,:) = sp_indices-twindow_si(1)-1;
spikewindow = cumsum(spikewindow,1);

spikewaves = hpData(spikewindow);
t_wave = [-twindow_si(1):twindow_si(2)]*si;

%%Plot filtered and unfiltered SOM recording
if nargin >= 3
    savename = varargin{1};
figure
    hold on
    plot([t_wave(1) t_wave(end)],[th th],'r--')
    plot([0 pt], [max(spikewaves(:))+2 max(spikewaves(:))+2])
    plot(t_wave,spikewaves,'Color',[0.5 0.5 0.5])
    text(0,max(spikewaves(:))+4,...
        [num2str(pt*10^3,3), ' +/- ',...
        num2str(pt_sd*10^3,3),' ms'])
    
    saveas(gcf,savename,'pdf')
end

end

