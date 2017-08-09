function [ viewinfo ] = EventVewPlot
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 

%Two options of current timepoint: Event Number, arbitraty timepoint
switch FO.viewmode
    case 'timepoint'
        timepoint = FO.currevent;
    case 'event'
        timepoint = FO.EventTimes(FO.currevent);
end


thiseventwin = timepoint+FO.winsize.*[-0.5 0.5];
inwinevents = FO.EventTimes(FO.EventTimes>=thiseventwin(1) &FO.EventTimes<=thiseventwin(2));

%Plot
%viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%set(gca,'ButtonDownFcn', @MouseClick)
hold(FO.viewwin,'off')
bz_MultiLFPPlot(FO.data.lfp,'timewin',thiseventwin,'spikes',FO.data.spikes,...
    'axhandle',FO.viewwin,'scaleLFP',FO.scaleLFP)
hold on
plot(FO.viewwin,inwinevents,zeros(size(inwinevents)),'o','color',[0 0.6 0])


%Passthrough info from the plot
viewinfo.inwinevents = inwinevents;
viewinfo.thiseventwin = thiseventwin;

end

