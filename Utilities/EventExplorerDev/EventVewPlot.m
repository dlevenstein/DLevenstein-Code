function [ viewinfo ] = EventVewPlot( FO,timepoint )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('timepoint','var')
    timepoint = FO.EventTimes(FO.currevent);
end

thiseventwin = timepoint+FO.winsize.*[-0.5 0.5];
inwinevents = FO.EventTimes(FO.EventTimes>=thiseventwin(1) &FO.EventTimes<=thiseventwin(2));

%Plot
%viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%set(gca,'ButtonDownFcn', @MouseClick)
bz_MultiLFPPlot(FO.data.lfp,'timewin',thiseventwin,'spikes',FO.data.spikes,'axhandle',FO.viewwin)
hold on
plot(FO.viewwin,inwinevents,zeros(size(inwinevents)),'o','color',[0 0.6 0])
hold off

%Passthrough info from the plot
viewinfo.inwinevents = inwinevents;
viewinfo.thiseventwin = thiseventwin;

end

