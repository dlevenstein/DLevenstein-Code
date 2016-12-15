function [ poiss ] = Poisson( r,s,dt )
%[ poiss ] = Poisson( r,s,dt ) returns the probability of observing s
%events in time window dt, given a poisson process or rate r.

poiss = ((r.*dt).^s).*exp(-r.*dt)./factorial(s);

end

