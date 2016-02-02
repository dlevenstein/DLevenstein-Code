function [trajectories] = GetTrajectories(spiketimes,basis,W,T,width,dt,varargin)
%GetTrajectories(spiketimes) projects the spikerate trajectories from 
%trials in spiketimes onto the same basis.
%
%Inputs
%   spiketimes: (N_trials x N_neurons) cell array, entry i,j is spiketimes 
%               of neuron j on trial i to match output format from 
%               TSObjects.
%   basis:      (N_neurons x N_dim) an orthonormal basis set - usually PCs
%   W:          [Wpre Wpost], added time on front and end of each trial
%   T:          (N_trials x 2) length vector with tstart and tend for each
%               trial
%   dt:         output time step
%   'plot'      (Optional) plots trajectories in top 3 dimensions
%   colorgroups (Optional - requires 'plot') assigns trajectories to groups
%               of same color
%
%Output structure 
%   trajectory(trial).traj
%   trajectory(trial).t
%   trajectory(trial).residual
%   trajectory(trial).group
%
%
%Dependencies:
%   SpktToRate(),SpktToSpkmat(),FConv(),Gauss()
%
%
%TO DO:
%   -Add Residual (amount of variance not explained by the projection)
%   -Bug: if projecting onto less than 3 dimensions... plot breaks
%   -Add group designation to output
%
%Last Updated: 5/20/15
%DLevenstein
%%
numtrials = length(spiketimes(:,1));
numneurons = length(spiketimes(1,:));
PLOTFIG = false;
%Plotting stuff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin)>0 
    assert(strcmp(varargin(1),'plot'),'First varargin must be "plot"')
	PLOTFIG = true;
    if length(varargin)>1
        %color by group
        colors = varargin{2};
        numgroups = length(unique(colors)); 
        colorpallet = 0.9*rand([numgroups,3]);
        colors = colorpallet(colors,:);
    else
        %each traj random color
        colors = 0.9*rand([numtrials,3]);       
    end
end

if PLOTFIG
    figure
    hold on
end

%Trajectories%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tr = 1:numtrials    
    %Convert time parameters to format for SpktToRate()
    T_trial = [W(1) T(tr,1) T(tr,2)-T(tr,1)+W(2)]; 
    %Get Trajectories    
    [spikerate,t] = SpktToRate(spiketimes(tr,:),width,T_trial,dt);
    trajectories(tr).traj = (basis'*spikerate')';
    trajectories(tr).t = t;    

%More Plotting Stuff%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if PLOTFIG 
        plot3(trajectories(tr).traj(1,1),trajectories(tr).traj(1,2),...
            trajectories(tr).traj(1,3),'ok',...
            trajectories(tr).traj(:,1),trajectories(tr).traj(:,2),...
            trajectories(tr).traj(:,3),'Color',colors(tr,:))
    end
end

if PLOTFIG
    grid on
    rotate3d
    xlabel('D1');ylabel('D2');zlabel('D3');
    hold off
end


end

