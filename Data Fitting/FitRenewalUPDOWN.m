function [ parmsHat,parmsCI,nloglik ] = FitRenewalUPDOWN( dwelltimes,SHOWFIG )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
if ~exist('SHOWFIG','var')
    SHOWFIG=false;
end


initalparms = [0 0.5 min(dwelltimes)];
lowerbound = [0 0 0];
upperbound = [1 1 (max(dwelltimes))];
opts = statset('MaxIter',50000, 'MaxFunEvals',100000,'FunValCheck','off');

[parmsHat,parmsCI] = mle(dwelltimes,'pdf',@RenewalFunction,...
    'options',opts,...%'optimfun','fmincon',...
    'start',initalparms,...
    'lowerbound',lowerbound,'upperbound',upperbound);

%acov = mlecov(parmsHat,dwelltimes,'pdf',@RenewalFunction);
%stderr = sqrt(acov);

nloglik = -sum(log(RenewalFunction(dwelltimes,parmsHat(1),parmsHat(2),parmsHat(3))));

%%
if SHOWFIG
    figure
    x = linspace(0,max(dwelltimes),100);
    [dwellhist,histbins] = hist(dwelltimes,20);
    dx = diff(histbins([1 2]));
    bar(histbins,dwellhist./(sum(dwellhist).*dx),'k')
    hold on
    plot(x,RenewalFunction(x,parmsHat(1),parmsHat(2),parmsHat(3)))
end

end

function pdf = RenewalFunction(data,p0,pinf,tau)
    p = (p0-pinf).*exp(-data./tau)+pinf;
    S = exp(-(p0-pinf).*tau.*(1-exp(-data./tau))-pinf.*data);
    
    pdf = S.*p;
end