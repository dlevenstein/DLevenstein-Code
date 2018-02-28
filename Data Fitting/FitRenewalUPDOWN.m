function [ parmsHat,parmsCI,nloglik ] = FitRenewalUPDOWN( dwelltimes,SHOWFIG )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
if ~exist('SHOWFIG','var')
    SHOWFIG=false;
end

if isnan(dwelltimes)
    parmsHat = [nan nan]; parmsCI = [nan nan]; nloglik= nan;
    return
end

%THREE PARM - p0,pinf,tau
% initalparms = [0 0.5 min(dwelltimes)];
% lowerbound = [0 0 0];
% upperbound = [1 1 (max(dwelltimes))];

%TWO PARM - pinf,tau
initalparms = [0.01 min(dwelltimes)];
lowerbound = [1e-6 1e-6];
upperbound = [1 10*(max(dwelltimes))];


opts = statset('MaxIter',1e6, 'MaxFunEvals',1e8,'FunValCheck','off');

[parmsHat,parmsCI] = mle(dwelltimes,'pdf',@RenewalFunction,...
    'options',opts,...%'optimfun','fmincon',...
    'start',initalparms,...
    'lowerbound',lowerbound,'upperbound',upperbound);

%acov = mlecov(parmsHat,dwelltimes,'pdf',@RenewalFunction);
%stderr = sqrt(acov);

nloglik = -sum(log(RenewalFunction(dwelltimes,parmsHat(1),parmsHat(2))));

%%
if SHOWFIG
    figure
    x = linspace(0,max(dwelltimes),100);
    [dwellhist,histbins] = hist(dwelltimes,20);
    dx = diff(histbins([1 2]));
    bar(histbins,dwellhist./(sum(dwellhist).*dx),'k')
    hold on
    plot(x,RenewalFunction(x,parmsHat(1),parmsHat(2)))
end

end

function pdf = RenewalFunction(data,pinf,tau)
%THREE PARM - p0,pinf,tau
%     p = (p0-pinf).*exp(-data./tau)+pinf;
%     S = exp(-(p0-pinf).*tau.*(1-exp(-data./tau))-pinf.*data);

%TWO PARM - pinf,tau
%    p = (0-pinf).*exp(-data./tau)+pinf;
%    S = exp(-(0-pinf).*tau.*(1-exp(-data./tau))-pinf.*data);
        
%TWO PARM - sigmoid
    p = pinf.*(data.^2)./(tau.^2+data.^2);
    S = exp(-pinf.*(data-tau.*atan(data./tau)));


    pdf = S.*p;
end