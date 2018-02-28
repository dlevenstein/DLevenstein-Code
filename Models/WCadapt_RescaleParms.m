function [ I,w,b ] = WCadapt_RescaleParms( I_rs,w_rs,b_rs,k,I0,tau )
%[ I,w,b ] = WCadapt_RescaleParms( I_rs,w_rs,b_rs,k,I0,tau ) takes values for the
%desired rescaled parameters and returns the actual model parameters


b = 4.*b_rs./k;
w = w_rs + 4.*(1+1./tau);
I = I_rs + I0 - 0.5.*(w-b);

end

