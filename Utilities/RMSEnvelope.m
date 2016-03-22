function [ R_data ] = RMSEnvelope( data, width, si )
%RMSEnvelope( data, width, si )
%
%TO DO:
%   -make work for vector of either orientation...
%
%Last Updated: 3/23/15

%Gaussian window
wind_x = [-2.5*width:si:2.5*width];
window = Gauss(wind_x,0,width);
window = window/sum(window);


%Squared
S_data = data.^2;

%Mean
%M_data = conv(S_data,window,'same');
M_data = FConv(window,S_data);

%Root
R_data = sqrt(M_data);

end

