function [  ] = ImportWMData(abffilename,savematfilename,cellnum,recnum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[raw_data,si,file_info]=abfload(abffilename);
SOM = raw_data(:,1);
LFP = raw_data(:,2);
EMG = raw_data(:,3);

%Extract relevant time information
si = si/(10^6);             %Convert sampling interval to s
sampfreq = 1/si;          %Sampling frequency
t_si = [0:length(LFP)-1];   %Time vector in samples
t = t_si*si;              %Time vector in seconds

%must be more than 3 minutes
%         if t(end) <= 60*3
%             keyboard
%         end


%Filter the EMG        
EMGrange = [400 3000];
EMG = FiltNPhase(EMG,EMGrange,sampfreq);

%Downsample for plots and LFP analysis frequency stuffs
downsample_factor = 20;     %target si: 1ms, sf: 1kHz
si_down = si*downsample_factor;
sf_down = 1/si_down;
SOM_down = downsample(SOM,downsample_factor);
LFP_down = downsample(LFP,downsample_factor);
EMG_down = downsample(EMG,downsample_factor);
t_down = downsample(t,downsample_factor);

%Isolate SOM Spikes from Cell Attached Recording 
%NOTE: Manually Improved by William, given better spikes in database
% savename = ['SpikeFigs/',...
%     'SOM',num2str(cellnum),'_',num2str(recnum),'_spikefig'];
% [SOM_spindices, SOM_spikes,t_wave,spikewaves,pt] = IsolateSpikes(SOM,...
%     sampfreq,savename);
% SOM_spindices_down = round(SOM_spindices/downsample_factor);
% SOM_spikes_down = zeros(size(SOM_down));
% SOM_spikes_down(SOM_spindices_down) = 1;

%Extract MUA from the LFP
[MUA_spindices,MUA_spikes,MUAPower,MUA,MUA_srate] = LFP2MUA(LFP,sampfreq);
MUA_spindices_down = round(MUA_spindices/downsample_factor);
%Take out overhanging indices and make a spike train vector
MUA_spindices_down(MUA_spindices_down<1|MUA_spindices_down>length(LFP_down))=[];
MUA_spikes_down = zeros(size(LFP_down));
MUA_spikes_down(MUA_spindices_down) = 1;

% LFPfilt = [0.1 300];
% LFP_down = FiltNPhase(LFP_down,LFPfilt,sf_down);

%Save to data structure for future speed of use.
recdata.t = t_down;
recdata.si = si_down;
recdata.si_raw = si;
recdata.LFP = LFP_down;
recdata.MUA = MUA_spindices*si;
recdata.EMG = EMG_down;
%recdata.spikes = SOM_spindices*si;
%recdata.peak_to_trough = pt;
recdata.cellnum = cellnum;
recdata.recnum = recnum;
    
    
%% Smooth the EMG and Identifying Whisking Epochs

%Parameters
gausswidth = 0.05;   %gaussian width for smoothing (s)
threshold = 0.7;    %EMG Threshold for Whisking (STDs)
minwhisk = 0.1;     %Minimum whisking length (s)
minNWh = 0.1;      %Minimum nonwhisking length (s)


%Z-Score the EMG and get EMG Envelope with RMS
EMGz = (EMG_down - mean(EMG_down))/std(EMG_down);
EMGsm = RMSEnvelope(EMGz',gausswidth,si_down);
%set EMG envelope minimum to 0
EMGsm = EMGsm-min(EMGsm);

%Smoothed EMG Stats
[EMGcounts,EMGvals] = hist(EMGsm,25);
EMGcounts = EMGcounts*si_down;

%Identify whisking on/offsets: smoothed EMG crosses threshold
wh_thresh = EMGsm > threshold;
wh_on = wh_thresh(2:end) > wh_thresh(1:end-1); 
wh_off = wh_thresh(2:end) < wh_thresh(1:end-1);
wh_on = find(wh_on == 1)+1;
wh_off = find(wh_off == 1)+1;

%If data starts/ends in the middle of an epoch, drop first/last trigger
if wh_off(1) < wh_on(1)
    wh_on = wh_on(1:end); wh_off = wh_off(2:end);
end
if wh_off(end) < wh_on(end)
    wh_on = wh_on(1:end-1); wh_off = wh_off(1:end);
end

%Drop nonwhisks epochs smaller than a minimum
[nwh_on,nwh_off] = MinEpochLength(wh_off,wh_on,minNWh,si_down);
%Remove whisking epochs smaller than minimum
[wh_on,wh_off] = MinEpochLength(nwh_off,nwh_on,minwhisk,si_down);
wh_on_t = wh_on*si_down;
wh_off_t = wh_off*si_down;
    

recdata.wh_on = wh_on_t;
recdata.wh_off = wh_off_t;
% recdata.srate_mean = srate_mean;
% recdata.srate_std = srate_std;

save(savematfilename,'recdata')
end

