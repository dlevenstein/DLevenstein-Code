function [ spikerate, t, spindices ] = SpktToRate(spiketimes,width,T,dt)
%SpktToRate(spiketimes,width,T,dt) Converts spike times to a continuous
%(t/dt x N_neurons) spikerate matrix.
%Note: spiketimes must be in seconds and dt must be multiple of 0.001
%
%
%%Inputs
%   spiketimes: (N_trials x N_neurons) cell array, entry i,j is spiketimes 
%               of neuron j on trial i to match output format from 
%               TSObjects.
%   T
%   T = [t_end]
%   T = [t_start t_end]              %For additional start time (t_start<0)
%   T = [t_start t_offset t_end]     %If t=0 is at aribtrary time t_offset
%
%
%
%Dependencies: 
%   FConv(),SpktToSpkmat(),Gauss()
%
%
%TO DO:
%   -Optional Arguement: smoothing kernal ('gauss','alpha','boxcar')
%   -Improve dt processing... Instead of going through SpktToSpkmat, add
%   gaussians centered on spiketimes?
%   -Can improve speed by doing "batch" convolution instead of looping each
%   cell through the fourier domain
%   -norm(spikerate(:,1)) should equal number of spikes... it does not.
%   -something is weird with t times after downsampling - found for 100ms
%   width with 10ms bin
%
%
%Last Updated: 5/18/15
%DLevenstein
%%
numneurons = length(spiketimes);

%Convert spiketimes to a spike matrix 
sm_dt = 0.001;  %dt for spike matrix
[spikemat,t,spindices] = SpktToSpkmat(spiketimes,T,sm_dt);


%The Gaussian kernel for smoothing
bound = ceil(4*width/sm_dt);    
wind_x = [-bound*sm_dt:sm_dt:bound*sm_dt]';
kernel = Gauss(wind_x,0,width);


%Smooth each neuron with the kernel
spikerate = zeros(size(spikemat));
for n = 1:numneurons
    spikerate(:,n) = FConv(kernel,spikemat(:,n));
end
%Downsample to desired dt.... should do this in a better way
downsamplefactor = dt/sm_dt;
spikerate = downsample(spikerate,downsamplefactor);
t = downsample(t,downsamplefactor);


end

