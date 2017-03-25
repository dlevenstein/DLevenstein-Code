function [ T,R,A,Inoise ] = WCAdapt_frange_run( simtime,dt,parms )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%PARMS (struct)
%   (Adaptation Parms)
%   .taus
%   .thresh
%   .gain
%   .mag
%   (Pop/Input Parms)
%   .W
%   .tau_r
%   .I_in
%   .noiseamp (white noise...)
%   .noisefreq
%
%
%Spring 2017 DLevenstein
%% PARMS
SHOWFIG = false;

taus = parms.taus;
thresh = parms.thresh;
gain = parms.gain;
mag = parms.mag;

W = parms.W;
I_in = parms.I_in;
tau_r = parms.tau_r;
noiseamp = parms.noiseamp;

ntaus = length(taus);
%% Initial Conditions: no activity

r_init = 0;
a_init = zeros(1,ntaus);

y0 = [r_init,a_init];       %combine initial conditions 
                            %into one long vector for ode45
                            
%time and solution arrays for ode45
t_tot = simtime/dt;             %number of iterations
tspan = [0:dt:simtime]';         %interval of integration

%% Noise for Input
Inoise = noiseamp*randn(length(tspan),1);
%Inoise = smooth(Inoise,1/(noisefreq*dt));
%http://www.mathworks.com/help/matlab/ref/ode45.html
%plot(tspan,Inoise)



%% System of Equations

    function dy = WCadapt_eqs(t, y)
%         if mod(t/simtime,0.05) == 0
%             display(num2str(t/simtime))
%         end
        %% separate the input vector into its e, i, and theta components
        r = y(1);
        a = y(2:end);
        
        %% THE DIFF.EQS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I_tot = W*r - sum(a) + I_in + interp1(tspan,Inoise,t)';
        
        F_I = 1./(1+exp(-I_tot));
        Ainf = mag./(1+exp(-(r-thresh)./gain));
        
        dr = (-r + F_I)./tau_r;
        da = (-a + Ainf')./taus';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% combine the outputs into one vector
        dy = [dr;da];
    end

%% RUN THE ODE SOLVER

[T,Y_sol] = ode45(@WCadapt_eqs,tspan,y0);

R = Y_sol(:,1);
A = Y_sol(:,2:end);


%% Plot the adaptation activation functions
if SHOWFIG

numr = 50;
r = linspace(0,1,numr);
for rr = 1:numr
    Ainf_R(:,rr) = mag./(1+exp(-(r(rr)-thresh)./gain));
end;


%pcolormap = [makeColorMap([1 1 1],[0.8 0.8 1],[0 0.5 0]);makeColorMap([0 0.5 0],[0.7 0 0],[0.8 0.5 0])];
adaptmap = [makeColorMap([0 0 0],[0 0 0.6],[0 0.6 0]);...
    makeColorMap([0 0.6 0],[0.6 0 0],[0.8 0.5 0.3])];

figure
colormap(adaptmap)
subplot(4,2,1)
    imagesc(log10(taus),r,Ainf_R')
    hold on
    plot(log10(taus),thresh,'k','linewidth',1)
    title(['AdaptMax: ',num2str(smax)])
    colorbar
    axis xy
    LogScale('x',10)

    Itot = linspace(-3,3,50);
subplot(4,2,2)
    plot(I_in,W,'*')
    hold on
    plot(Itot,1./(1+exp(-Itot)),'k','linewidth',1);
    xlim([-3 3]);ylim([0 8])
    
subplot(4,1,2)
    plot(T,R,'k','linewidth',1)
    axis tight
subplot(4,1,3)
    imagesc(T,log10(taus),A')
    axis xy
    LogScale('y',10)
    %caxis([0 1])
    colorbar('east')
subplot(4,1,4)
    plot(T,Inoise,'color',0.5.*[1 1 1])
    hold on
    plot(T,-sum(A,2),'linewidth',2)
   % hold on
    plot(T,W.*R,'linewidth',2)

    axis tight

 %   figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/SlowOscillation/Modeling/Figures';
%NiceSave('simdata',figfolder,'WCAdapt_frange')

end

end

