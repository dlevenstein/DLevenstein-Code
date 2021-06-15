function [ phandle ] = StateScorePlot( stateints,colors,varargin )
%UStateScorePlot( stateints,colors ) plots a state score indicator plot
%from state intervals with designated colors.
%
%INPUTS
%   stateints   structure with intervals or 
%               {Nstates} cell array of [Nints x 2] time intervals
%   colors      {Nstates} cell array of colors
%%
% parse args
p = inputParser;
addParameter(p,'y','top')
addParameter(p,'LineWidth',8)
parse(p,varargin{:})
yplot = p.Results.y;
LineWidth = p.Results.LineWidth;
%%
plotylimits = get(gca,'ylim');
if strcmp(yplot,'top')
    yplot = plotylimits(2)*1.05;
    ytop = yplot;
else
    ytop = plotylimits(2);
end
yrange = plotylimits(2)-plotylimits(1);


if isstruct(stateints)
    stateints = struct2cell(stateints);
end

numstates = length(stateints);
for ss = 1:numstates
    if isa(stateints{ss},'intervalSet')
        stateints{ss} = [Start(stateints{ss},'s'), End(stateints{ss},'s')];
    end
end

yscale = {0.11,0.06,0.01};
%statey = cellfun(@(X) yplot+yrange*X,yscale,'UniformOutput',false);


%%
hold on
    for ss = 1:numstates
        plot(stateints{ss}',yplot*ones(size(stateints{ss}))','Color',colors{ss},'LineWidth',LineWidth)
        %alpha(0.5)
    end

ylim([plotylimits(1) ytop]);  
end

