function [ output_args ] = StateScorePlot( stateints,colors )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
yrange = get(gca,'ylim');
ylow = yrange(2)*1.01;
yrange = yrange(2)-yrange(1);

numstates = length(stateints);
for ss = 1:numstates
    if isa(stateints{ss},'intervalSet')
        stateints{ss} = [Start(stateints{ss},'s'), End(stateints{ss},'s')];
    end
end

yscale = {0.11,0.06,0.01};
statey = cellfun(@(X) ylow+yrange*X,yscale,'UniformOutput',false);


%%
    for ss = 1:numstates
        plot(stateints{ss}',ylow*ones(size(stateints{ss}))',colors{ss},'LineWidth',10)
    end

end

