function [ phandle ] = StateScorePlot( stateints,colors )
%UStateScorePlot( stateints,colors ) plots a state score indicator plot
%from state intervals with designated colors.
%
%INPUTS
%   stateints   {Nstates} cell array of [Nints x 2] time intervals
%   colors      {Nstates} cell array of colors
%%
plotylimits = get(gca,'ylim');
ylow = plotylimits(2)*1.05;
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
statey = cellfun(@(X) ylow+yrange*X,yscale,'UniformOutput',false);


%%
hold on
    for ss = 1:numstates
        plot(stateints{ss}',ylow*ones(size(stateints{ss}))','Color',colors{ss},'LineWidth',8)
    end

ylim([plotylimits(1) ylow]);  
end

