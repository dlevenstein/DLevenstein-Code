function [ ccghist,ccgbins ] = EventCCG(events1,events2,timeparms)
%EventCCG(events1,events2)
%
%INPUTS
%   events1/2   vector of event times or [Nints x 2] interval start and
%               end times.
%   timeparms.
%       twin      (s) (before/after)
%       numbins   (before/after)
%       normbins  (within start-end)
%
%OUPUTS
%               if event times, calculates standard event CCG (times of
%               event 2 relative to event 1).
%               if events 1 or 2 are intervals - calculates internal
%               time-normalized CCG where start/end times are calculated
%               relative to start/end times and internally time normalized
%               time (i.e. a time axis looks like:
%               -1s---S----E---+1
%
%
%TO DO
%-inputParser
%
%DLevenstein 2016
%% Parms (to be made input options with inputParser)
twin = timeparms.twin; 
numbins = timeparms.numbins; %
binsize = (twin*2)./numbins; 
normbins = timeparms.normbins;  %(within start-end)
%NOTE: need to renormalize histogram to counts/time


%%
if isrow(events1) && length(events1)>2
    events1 = events1';
end
if isrow(events2) && length(events2)>2
    events2 = events2';
end

if isempty(events1)
    ccgbins = linspace(-twin,twin,numbins);
    ccghist = nan(size(ccgbins));
%    keyboard
    return
end

event2type = length(events2(1,:));
event1type = length(events1(1,:));

numE1 = length(events1(:,1));
numE2 = length(events2(:,1));

%%

switch event2type
    case 1
        
        %The meat of calculating CCG is here: 
        %Calculates relative time to start and end of events1.
        relativetime = [];
        for ee = 1:numE1
            relevantwindow = events1(ee,:) + twin*[-1 1];
            relevantevents = events2>relevantwindow(1) & events2<relevantwindow(2);
            relativetime = [relativetime; bsxfun(@(X,Y) X-Y,events2(relevantevents),events1(ee,:))];
        end
        
        switch event1type
            
            %Standard CCG
            case 1 
                ccgbins = linspace(-twin,twin,numbins+1);
                ccgbins = ccgbins(1:end-1)+0.5.*diff(ccgbins([1 2]));
                [ccghist] = hist(relativetime,ccgbins);
                ccghist = ccghist./(binsize*numE1);
                
                figure
                    bar(ccgbins,ccghist)
                    xlabel('t lag (s)');ylabel('Event Rate (Events/s)')
                
            %If event1 has start and end times        
            case 2 
                %Identify internal events2 and calculate relative time.
                internalevents = sign(relativetime(:,1))==1 & sign(relativetime(:,2))==-1;
                internaltime = relativetime(internalevents,1)./(relativetime(internalevents,1)-relativetime(internalevents,2));
                %Histogram of internal events.
                internalbins = linspace(0,1,normbins+1);
                internalbins = internalbins(1:end-1)+0.5.*diff(internalbins([1 2]));
                [internalhist] = hist(internaltime,internalbins);
                totalinternaltime = sum(events1(:,2)-events1(:,1));
                internalhist = internalhist./(totalinternaltime./normbins);

                %External Events
                relativetime = relativetime(~internalevents,:);
                [~,whichone] = min(abs(relativetime),[],2);
                whichone = sub2ind(size(relativetime),1:length(whichone),whichone');
                relativetime = relativetime(whichone);
                
                ccgbins = linspace(-twin,twin,numbins+1);
                ccgbins = ccgbins(1:end-1)+0.5.*diff(ccgbins([1 2]));
                [ccghist] = hist(relativetime,ccgbins);
                ccghist = ccghist./(binsize*numE1);
                ccgbins(ccgbins>0) = ccgbins(ccgbins>0)+1;
                
                ccgbins = [ccgbins,internalbins];
                ccghist = [ccghist,internalhist];
                
                figure
                    bar(ccgbins,ccghist)
                    hold on
                    plot([0 0],get(gca,'ylim'),'k')
                    plot([1 1],get(gca,'ylim'),'k')
                    xlabel('t lag (s)');ylabel('Event Rate (Events/s)')
                    set(gca,'XTick',[-twin:twin+1])
                    set(gca,'XTickLabels',{-twin:-1,'S','E',1:twin})
        end
        

    %Calculates CCG for event2 start and end times separately
    case 2  

        [event2startccg,ccgbins] = EventCCG(events1,events2(:,1));
        [event2endccg,ccgbins] = EventCCG(events1,events2(:,2));
        ccghist = [event2startccg;event2endccg];

        figure    
            bar(ccgbins,event2startccg,'g')
            hold on
            bar(ccgbins,-event2endccg,'r')
            plot(ccgbins,event2startccg-event2endccg,'k','LineWidth',1)
            plot([0 0],max(abs(get(gca,'ylim')))*[-1 1],'k')
            xlabel('t lag (s)');ylabel({'Event Rate (Events/s)','Ends <----------> Starts'})
end




end

