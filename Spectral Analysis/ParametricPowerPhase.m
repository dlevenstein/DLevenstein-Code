function [ R ] = ParametricPowerPhase(rho,theta,alpha,beta,phi,pthresh,thmod)
%[R] = ParametricPowerPhase(rhovals,thetavals,alpha,beta,phi,pthresh,thmod)
%Parametric function for power-modulated phase coupling.
%
%PARAMETERS
%   alpha       rate-power modulation
%   beta        baseline rate
%   phi         preferred phase 
%   pthresh     phase-coupling power threshold
%   thmod       phase modulation strength
%% PowerPhaseCoupling

% numbins = 50;
% rhovals = linspace(0,1,50)';
% thetavals = linspace(0,2.*pi,50);
% 
% alpha = 0.2;
% beta = 1;
% phi = 0.2;
% 
% 
% pthresh = 0.5;
% thmod = 1;

%%

%rho = repmat(rhovals,1,numbins);
%theta = repmat(thetavals,numbins,1);

%Phase Modulation Magnitude
A = thmod.*(rho-pthresh)./(1-pthresh);
A(rho<pthresh) = 0;

%Parametric Power/Phase Modulation Function
R = alpha.*rho+beta+A.*sin(theta+phi);

%Rectification
R(R<0)=0;

%%
% figure
% imagesc(thetavals,rhovals,R)
% hold on
% imagesc(thetavals+2.*pi,rhovals,R)
% axis tight
% axis xy
% colorbar
% xlabel('Phase');ylabel('Power')



end

