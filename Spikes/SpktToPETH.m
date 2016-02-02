function [PETHspikes,PETHrates,t] = SpktToPETH(spiketimes,T,gausswidth,eventsort,cellsort,dt)
%SpktToPETH Summary of this function goes here
%
%Inputs
%   spiketimes  (N_events x N_neurons) cell array, entry i,j is spiketimes 
%               of neuron j on trial i to match output format from 
%               TSObjects.
%   T:          (N_events x 2) length vector with tstart and tend for each
%               trial.
%   gausswidth  width of gaussian window for firing rate smoothing
%   eventsort   order to index events for PETHspikes
%               example: [events_sorted,eventsort] = sort(eventproperty);
%   cellsort    order to index cells for PETHspikes
%   dt          time step for rate vector output
%
%
%Outputs
%   PETHspikes  (N_spikes x 2) matrix for plotting spikes 
%               PETHspikes(1,:) is spiketimes; PETHspikes(2,:) is 
%               event/cell index
%   PETHrates   (t x N_cells) matrix, column i is the event-averaged 
%               firing rate of cell i
%   t           event time vector (s)
%
%
%TO DO:
%   -Add divider lines for plot output - see UPOnset.m
%   -Change PETHspikes to PETHspikeplot
%       PETHspikeplot.spiketimes
%       PETHspikeplot.spikeindices
%       PETHspikeplot.divlines
%       etc...
%   -Add option to not need sort for events?
%   -BUG: eventsort must be column vector
%
%Last Updated: 5/20/15
%DLevenstein

%%
numtrials = length(spiketimes(:,1));
numneurons = length(spiketimes(1,:));
Tlength = T(:,2)-T(:,1);

%Align all spike times to event start
Tstarts_cell = num2cell(repmat(T(:,1),1,numneurons));
spiketimes_aligned = cellfun(@minus,spiketimes,Tstarts_cell,'UniformOutput',false);
PETHspikes = vertcat(spiketimes_aligned{:}); %Spiketimes output vector

%Collapse spikes from each cell to a single spiketimes cell array "event"
spiketimes_allevents = cell(1,numneurons);
for c = 1:numneurons
    spiketimes_allevents{c} = vertcat(spiketimes_aligned{:,c});
end

%Calculate Spike Rate for each cell
T_rate = [0 max(Tlength)];
[PETHrates,t] = SpktToRate(spiketimes_allevents,gausswidth,T_rate,dt);

%Divide all-event rate by the number of events lasting up to each time 
%point... get rid of this for loop, do this better 
eventvec = zeros(size(t));
for d = 1:length(t)
    eventvec(d) = sum(Tlength>t(d));
end
PETHrates = PETHrates./repmat(eventvec,[1,numneurons]);

[dum,eventindex] = sort(eventsort);
[dum,cellindex] = sort(cellsort);
%Spike Index for PETHspikes sorted by event and cell
eventindex = num2cell(repmat(eventindex,1,numneurons)); %By UP state length
cellindex = num2cell(repmat(cellindex-1,numtrials,1).*numtrials);  %by firing rate
spikeindex = cellfun(@(x) x.*0,spiketimes_aligned,'UniformOutput',false);
spikeindex = cellfun(@plus,spikeindex,eventindex,'UniformOutput',false);
spikeindex = cellfun(@plus,spikeindex,cellindex,'UniformOutput',false);
PETHspikes = [PETHspikes,vertcat(spikeindex{:})];

end

