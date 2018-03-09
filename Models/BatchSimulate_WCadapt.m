function [ dwelltimes,ratehist,fspec ] = BatchSimulate_WCadapt(modelparms,simparms)
%[dwelltimes,ratehist] = BatchSimulate_WCadapt(modelparms,simparms)
%simulates the adapting wilson-cowan model for a number of parameter
%values in parallel.
%
%INPUT
%   modelparms(n)   is a structure with the parameters for each of the n
%                   simulations
%       .N_neurons
%       .I_in
%       .W
%       .beta
%       .tau_r
%       .tau_a
%       .A0
%       .Ak
%       .noiseamp
%       .noisefreq
%   simparms(n)
%       .simtime
%       .dt
%
%
%%

numsims = length(modelparms);
timespent = 0;

parms = modelparms;

simtime = simparms.simtime;
dt = simparms.dt;

dwelltimes = struct('UP',[],'DOWN',[]);
dwelltimes(numsims)=dwelltimes;
ratehist = struct('bins',[],'hist',[],'mean',[],'std',[],'thresh',[]);
ratehist(numsims) = ratehist;

numfreqs = 200;
fspec = struct('freqlist',[],'spec',[]);
fspec(numsims) = fspec;

parfor_progress(numsims);
tstart = tic;
%% Loop Simulation
parfor nn = 1:numsims

    tempstruct_rate = struct('bins',[],'hist',[],'mean',[],'std',[],'thresh',[]); %temporary structure to hold rate distribution
    tempstruct_dwell = struct('UP',[],'DOWN',[]); %temporary structure to hold dwelltime distribution
    
    timespent=toc(tstart);
    percdone = parfor_progress;

    estimatedtotal = timespent./(percdone./100);
    estimatedremaining = estimatedtotal-timespent;
    display(['Percent Done: ',num2str(percdone),...
        '.  Time Spent: ',num2str(round(timespent./60,1)),...
        '.  Est. Total Time: ',num2str(round(estimatedtotal./60,1)),...
        'min.  ETR: ',num2str(round(estimatedremaining./60,1)),'min.'])
    
    %% Run the Simulation

    [ ~, Y_sol ] = WCadapt_run(simtime,dt,parms(nn));
    r = Y_sol(:,1);
    %a = Y_sol(:,2);

    %% Calculate Rate distribution
    numratebins = 50;
    tempstruct_rate.bins = linspace(0,1,numratebins);
    tempstruct_rate.hist = hist(r,tempstruct_rate.bins);
    tempstruct_rate.mean = mean(r);
    tempstruct_rate.std = std(r);
    
    %% Calculate DWELL times
    [tempstruct_rate.thresh,cross,~] = BimodalThresh(r,'Schmidt');
    
    if isempty(cross.upints) || length(cross.upints) <=2
        tempstruct_dwell.UP = nan;
        tempstruct_dwell.DOWN = nan;
    else
        tempstruct_dwell.UP = cross.upints(:,2)-cross.upints(:,1);
        tempstruct_dwell.DOWN = cross.downints(:,2)-cross.downints(:,1);
    end
    
    ratehist(nn) = tempstruct_rate;
    dwelltimes(nn) = tempstruct_dwell;

    %% Calculate Frequency Spectrum
    
%     fspec(nn).freqlist = logspace(-3,1,numfreqs);
%     window = 5000;
%     noverlap = 2500;
%     window = window/dt;
%     noverlap = noverlap/dt;
%     [FFTspec,FFTfreqs,t_FFT] = spectrogram(r,window,noverlap,fspec(nn).freqlist,1/dt);
%     fspec(nn).spec = mean(abs(FFTspec),2)';

end

parfor_progress(0);

end

