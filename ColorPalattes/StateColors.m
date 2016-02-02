function colors = StateColors
% Generates RGB Triplets based on red-blue continuum:
% 
% INPUT
% -numcolors: number of colors you want out
% OUTPUT
% - colors a 3 column matrix signifying Red Green and Blue values with each
% column, number of rows = numcolors
% 
% Brendon Watson 2015

colors(1,:) = [0 0 0];
colors(2,:) = [.3 .4 1];
colors(3,:) = [.8 0 0];

% r = linspace(.8,.3,numcolors)';
% g = linspace(0,.4,numcolors)';
% b = linspace(0,1,numcolors)';

% colors = [r g b];
