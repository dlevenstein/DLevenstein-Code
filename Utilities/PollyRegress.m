function [ betas_opt, error ] = PollyRegress( x, y, order )
%PollyRegress( x, y, order ) Uses  least squares regression to fit data y
%as a polynomial function of explanatory variable x. Returns the values
%of beta that minimize error and the error as [betas_opt, error].
%
%First made for Math Tools HW3
%Last Updated: 3/31/14
%DLevenstein

%Make the Explanatory Variable matrix for a polynomial regression of
%appropriate order 
%Note: there's definitely a more clever way to do this without a loop...
X_exp = x.^0;
    for i = 1:order
        X_exp = [X_exp,x.^i];
    end

%Use the SVD to solve the values of beta that minimize ||X*beta-y||^2.
%i.e. beta_opt = ((X'X)^-1)*X'*y
%[U,S,V] = svd(X_exp);
%betas_opt = inv(V*S'*S*V')*X_exp'*y;
%Note: the easy way to do this is
betas_opt = pinv(X_exp)*y;
%but the homework instructions say to use svd...

error = (X_exp*betas_opt-y)'*(X_exp*betas_opt-y);

end

