function [ whisk ] = GetWhiskFromEMG( baseName,basePath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% DEV
basePath = '/mnt/proraidDL/Database/WMProbeData/';
baseName = 'Layers_LFP_Test02_170323_151411';

%%
abfname = fullfile(basePath,baseName,[baseName,'.abf']);
analogName = fullfile(basePath,baseName,['analogin.dat']);

%%
%% Clampex File
timechan = 1;
emgchan = 2;
sf_abf = 20000; %Sampling Frequency of the .abf file


abffile = abfload(abfname);
pulse_abf = abffile(:,1);
EMG = abffile(:,2);

t_abf = [1:length(EMG)]'./sf_abf;

%%
% figure
% plot(t_abf,EMG,'k')
%% Get the whisking envelope
EMGrange = [400 3000];
EMG = FiltNPhase(EMG,EMGrange,sf_abf);

downsamplefactor = 16; %Downsample to same as the LFP;
EMG = downsample(EMG,downsamplefactor);
t_EMG = downsample(t_abf,downsamplefactor);
sf_down = sf_abf./downsamplefactor;
%%

EMGparms.gausswidth = 0.05;  %Gaussian width for smoothing (s)
EMGparms.threshold = 2.5;    %EMG Threshold for Whisking (modSTDs)
EMGparms.minwhisk = 0.1;     %Minimum whisking duration (s)
EMGparms.minNWh = 0.1;       %Minimum nonwhisking duration (s)

%% Z-Score the EMGZ and get EMG envelope with RMS
EMGz = NormToInt(EMG,[],[],'modZ'); %Modified Z score - robust to outliers
EMGsm = RMSEnvelope(EMGz,EMGparms.gausswidth,1/sf_down);
EMGsm = EMGsm-min(EMGsm);

whisk.EMGenvelope = EMGsm;

%% Identify Whisking on/offsets: EMG envelope crosses threshold
wh_thresh = EMGsm > EMGparms.threshold;
wh_on = find(wh_thresh(2:end)>wh_thresh(1:end-1))+1; %whisking onsets (si)
wh_off = find(wh_thresh(2:end)< wh_thresh(1:end-1))+1;%whisking offsets (si)

%% If data starts/ends in the middle of an epoch, drop first/last trigger
if wh_off(1)<wh_on(1)
    wh_off = wh_off(2:end);
end
if wh_off(end) < wh_on(end)
    wh_on = wh_on(1:end-1);
end

%Drop nonwhisk epochs smaller than a minimum
[nwh_on,nwh_off] = MinEpochLength(wh_off,wh_on,EMGparms.minNWh,1/sf_down);
%Drop whisking epochs smaller than a minimum
[wh_on,wh_off] = MinEpochLength(nwh_off,nwh_on,EMGparms.minwhisk,1/sf_down);
%%
figure
subplot(2,1,1)
plot(t_EMG,EMGz,'k')

hold on
plot(t_EMG,EMGsm,'b')

subplot(2,1,2)
plot(t_EMG,EMGmodz,'k')

hold on
plot(t_EMG,EMGsm_modz,'b')

%%
figure




%% Identify Pulses from Camera
pulsethreshold =0.5;  %Adjust this later to set based on input.
pulseonsets = find(diff(timepulses<pulsethreshold)==1);
pulset = t_pulse(pulseonsets);


pulset(1) = []; %remove the first trigger... make this more rigorous later 

end

