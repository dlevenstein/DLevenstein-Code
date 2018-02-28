function [ s ] = PoissonRateSpikeBins(r,dt,T)
%[ s ] = PoissonRateSpikeBins(r,dt,T)
%Given a (NNeurons) vector of rates and a time window dt, this function
%returns a probabilisitically determined spike matrix (1 = spike in time
%window) for the population.
%
%DLevenstein 2017
%% DEV
%dt = 0.001;
%r = [0.1 1 10 100]';
%T =1000;

%%
if iscolumn(r)
    r = r';
end

if exist('T','var')
    r = repmat(r,T,1);
end

%Probability of spiking in time window of size dt for each neuron
p_spike = r.*dt;

if any(p_spike>=1)
    display(['FYI: Rate of >=1 of your neurons is too high for your time bin...'...
        'automatically setting to 1 spike'])
end

%
s = rand(size(r));
s = s<p_spike;


end

