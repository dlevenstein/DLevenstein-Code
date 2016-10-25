function [lognorm] = LogNorm(x,mu,sig)
%LogNorm(x,mu,s) returns a normalized lognorm distirbution with location mu 
%and scale sig

p1 = -0.5.*((log(x)-mu)./sig).^2;   %Term in the exponent
p2 = (x.*sig .* sqrt(2.*pi));          %Normalization

lognorm = exp(p1) ./ p2;          %Put it all together!

end
