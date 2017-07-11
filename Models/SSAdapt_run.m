function [ T, Y_sol,Inoise ] = SSadapt_run(simtime,dt,parms)
%UNTITLED2 Summary of this function goes here
%
%parms
%   .N_neurons
%   .W
%   .I_in
%   .beta
%   .tau_r
%   .tau_a
%   .Ak
%   .A0
%   .noiseamp
%   .noisefreq
%   .pulseparms (optional)   [t_onset magnitude duration]
%
%
%TO DO:
%   -Allow for general activation function as input
%
%%
parms.samenoise = false;

%% Establish Vectors etc
N_neurons = parms.N_neurons;
W = parms.W;
I_in = parms.I_in;
beta = parms.beta;
tau_r = parms.tau_r;
tau_a = parms.tau_a;
Ak = parms.Ak;
A0 = parms.A0;
k = parms.k;
h = parms.h;
n = parms.n;
noiseamp = parms.noiseamp;
noisefreq = parms.noisefreq;

pulseparms = [0 0 0];
if isfield(parms,'pulseparms')
    if isempty(parms.pulseparms);pulseparms = [0 0 0];end
    pulseparms = parms.pulseparms;
end


%initial conditions: no activity. 
%                    -can add different initial conditions if desired.
%r_init = zeros(1,N_neurons);
r_init = rand(1,N_neurons);
%a_init = zeros(1,N_neurons);
a_init = rand(1,N_neurons);

y0 = [r_init,a_init];       %combine initial conditions 
                             %into one long vector for ode45
                                            
%time and solution arrays for ode45
t_tot = simtime/dt;             %number of iterations
tspan = [0:dt:simtime]';         %interval of integration


% %Do these need to be here?....
% r_sol = zeros(t_tot+1,N_neurons);       %solution array for E cells
% a_sol = zeros(t_tot+1,N_neurons);     %solution array for theta  
% T = tspan;                     %column vector of time points
% Y_sol = [r_sol,a_sol];          %solution array for ode45

%% Noise for Input
%noisefilter = [0.5 10];
%Inoise = noiseamp.*wgn(length(tspan),1,0);
switch parms.samenoise
    case true
        Inoise = noiseamp*randn(length(tspan),1);
        Inoise = repmat(Inoise,N_neurons,1);
    case false
        Inoise = noiseamp*randn(length(tspan)*N_neurons,1);
end
%Inoise = FiltNPhase(Inoise,noisefilter,1000);
Inoise = smooth(Inoise,(1/noisefreq)*dt);
Inoise = reshape(Inoise,[length(tspan),N_neurons]);
%http://www.mathworks.com/help/matlab/ref/ode45.html
%plot(tspan,Inoise)
%% System of Equations

    function dy = SSadapt_eqs(t, y)
%         if mod(t/simtime,0.05) == 0
%             display(num2str(t/simtime))
%         end
        %% indices
        r_indexlow = 1;
        r_indexhigh = N_neurons;
        a_indexlow = r_indexhigh + 1;
        a_indexhigh = r_indexhigh + N_neurons;
        
        %% separate the input vector into its e, i, and theta components
        r = y(r_indexlow:r_indexhigh);
        a = y(a_indexlow:a_indexhigh);
        
        %% THE DIFF.EQS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_tot = W*r - a + I_in + interp1(tspan,Inoise,t)' + pulsefun(t,pulseparms);
        
        F_I =  k.*(heaviside(I_tot-h).*(I_tot-h)).^n;
        Ainf = beta./(1+exp(-Ak.*(r-A0)));
        
        dr = (-r + F_I)./tau_r;
        da = (-a + Ainf)./tau_a;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% combine the outputs into one vector
        dy = [dr;da];
    end

%% RUN THE ODE SOLVER

[T,Y_sol] = ode45(@SSadapt_eqs,tspan,y0);


%% Function for pulse
function Ipulse_out = pulsefun(t_in,pulseparms_in)
    t_onset=pulseparms_in(:,1);magnitude=pulseparms_in(:,2);duration=pulseparms_in(:,3);
    t_offset = t_onset + duration;
    
    Ipulse_out = (t_in >= t_onset & t_in <= t_offset).*magnitude;
%     if t_in < t_onset
%         Ipulse_out = 0;
%     elseif t_in >= t_onset && t_in <= t_offset
%         Ipulse_out = magnitude;
%     elseif t_in > t_offset
%         Ipulse_out = 0;
%     end
end

end

