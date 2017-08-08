function [ EEoutput ] = EventExplorer(events,basePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%DLevenstein 2017
%%
%eventsname = 'SlowWave';
%basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);


%%

%A0 = EE_Initiate(
%% Select the Events
%eventsfile = fullfile(basePath,[baseName,'.',eventsname,'.events.mat']);

%events = bz_LoadEvents(eventsname);

if isstring(events)
    events = bz_LoadStates(basePath,eventsname);
    eventstype = 'states';

    switch eventstype
        case 'states'
            intnames = [fieldnames(events.ints);fieldnames(events)];
            [s,v] = listdlg('PromptString','Which interval type to explore?',...
                'SelectionMode','single','ListString',intnames);
            exploreintname = intnames{s};
            try
            exploreint = events.ints.(exploreintname);
            catch
                exploreint = events.(exploreintname);
            end
    end

elseif isstruct(events)
    exploreint = events.timestamps;
    exploreintname = events.detectorinfo.detectorname;
end
%% 
FO.baseName = baseName;
FO.EventTimes = exploreint;
FO.EventName = exploreintname;
FO.basePath = basePath;
%% Make the EventExplorer Figure
%EE_Initiate( FO )
usechans = events.detectorinfo.detectionparms.SWchannel;
[ SleepState ] = bz_LoadStates(FO.basePath,'SleepState');
restrictint = SleepState.ints.NREMstate;
%lfp = bz_GetLFP(usechans,'basepath',FO.basePath,'intervals',restrictint);
lfp = bz_GetLFP(usechans,'basepath',FO.basePath);
spikes = bz_GetSpikes('basepath',FO.basePath);
%%

posvar = get(0,'Screensize');
posvar(1) = 20;
posvar(2) = 20;
posvar(3) = posvar(3)-100;
posvar(4) = posvar(4)-100;

%Old figure should be closed
oldfig = findobj('tag','EventExplorerMaster');
close(oldfig)
%Start the figure
FO.fig = figure('Position', posvar);
set(FO.fig, 'numbertitle', 'off', 'name', ['Recording: ', FO.baseName,'. Events: ',FO.EventName]);
set(FO.fig, 'Tag', 'EventExplorerMaster');
%set(FO.fig, 'CloseRequestFcn', {@CloseDialog}); The close function
set(FO.fig,'WindowButtonDownFcn', {@MouseClick});

%From StateEditor
% FO.fig = figure('KeyReleaseFcn', {@DefKey}, 'Position', posvar);
% set(FO.fig, 'numbertitle', 'off', 'name', ['Event: ', FO.baseName]);
% set(FO.fig,'WindowButtonDownFcn', {@MouseClick}, 'WindowButtonUpFcn', {@unMouseClick}, 'Units', 'normalized');
% set(FO.fig, 'WindowButtonMotionFcn', {@Nothing}, 'WindowScrollWheelFcn', {@MouseScroll});
% set(FO.fig, 'CloseRequestFcn', {@CloseDialog});
% set(FO.fig, 'Tag', 'EventExplorerMaster');


%Store the data in the figure - do this at the end of each function?
guidata(FO.fig, FO);


%% Select the time windows to look at 
winsize = 8; %(s) user input this or determine based on event frequency
numevents = 25; %number of events to look at. determine to maximize sampling or have user input (FO.uinput.field?)
%Selecting from events:
randevents = randsample(FO.EventTimes,numevents);
%Selecting from random times
[NREMtime,~,~] = RestrictInts(lfp.timestamps,restrictint);
randevents = randsample(NREMtime,numevents);
%Find any events that are within winsize of another event to remove them?
%Also should look at windows in which no events were detected?  Maybe just
%random windows in the detectionwin
closeevents = abs(diff(randevents))<winsize;
% numtries =0;
% while any(closeevents) && numtries<50
%     randevents(closeevents)
% end
%Problem:overlapping ints
%[ ints ] = MergeSeparatedInts( ints,minseparation )

%% Loop figure to show events and get user input

%Store the current user action
FO.currentuseraction = 'MarkEvents';
FO.markedevents = [];
guidata(FO.fig, FO);

%Set up the counters
counter = numevents;
miss=[];correct=[];falsealarm=[];


correcttext = uicontrol(FO.fig,'Position',[40 60 125 10],'style','text',...
    'string','Correct: 0','HorizontalAlignment','left'); 
misstext = uicontrol(FO.fig,'Position',[40 45 125 10],'style','text',...
    'string','Miss: 0','HorizontalAlignment','left');
FAtext = uicontrol(FO.fig,'Position',[40 30 125 10],'style','text',...
    'string','FA: 0','HorizontalAlignment','left');
countetext = uicontrol(FO.fig,'Position',[40 80 125 10],'style','text',...
    'string',['Windows Remaining: ',num2str(counter)],'HorizontalAlignment','left');


QUITLOOP = false;
while counter>0 && QUITLOOP~=true
    %EventUI(FO)   %function that will run the user interface
    thiseventtime = randevents(counter);
    thiseventwin = thiseventtime+winsize.*[-0.5 0.5];
    inwinevents = FO.EventTimes(FO.EventTimes>=thiseventwin(1) &FO.EventTimes<=thiseventwin(2));
    
    lfpwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
    %set(gca,'ButtonDownFcn', @MouseClick)
    bz_MultiLFPPlot(lfp,'timewin',thiseventwin,'spikes',spikes)
    hold on
    plot(inwinevents,zeros(size(inwinevents)),'o','color',[0 0.6 0])
    
    %Get user feedback
    title('Left click missed events. Right click incorrectly ID''d events.')
    %[x,y,button] = ginput;

    %if any(button==27);QUITLOOP=true;end
    
    %Update window to show how many misses/FAs have been ID'd
    %When click.... mark the event as FA or miss on screen... with red o,x
    %Top window to give context?
    %Instructions window
    
    %Buttons for next loop or finishing
    h = uicontrol('Position',[380 20 125 40],'String','Next Window (Return)',...
              'Callback','uiresume(gcbf)');
    qutbtn = uicontrol('Position',[380 65 125 40],'String','Quit Early (Esc)',...
             'Callback','QUITLOOP=true;uiresume(gcbf)');
    
    uiwait(FO.fig)
    hold off
    counter = counter-1;
    
    %Tally the user selections for that window
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj);
    if ~isempty(FO.markedevents)
        miss = [miss; FO.markedevents(FO.markedevents(:,3)==1,1)];
        
        %This is messy to account for interp errors with 0-1 reference points
        if isempty(FO.markedevents(FO.markedevents(:,3)==3,1))
            newfalsealarms = [];
        elseif length(FO.markedevents(FO.markedevents(:,3)==3,1))==1 && length(inwinevents)==1
            newfalsealarms = inwinevents;
        else
        newfalsealarms = interp1(inwinevents,inwinevents,FO.markedevents(FO.markedevents(:,3)==3,1),'nearest');
        end
        falsealarm = [falsealarm; newfalsealarms];
        
    end
    correct = [correct; inwinevents(~ismember(inwinevents,falsealarm))];
    FO.markedevents = [];  guidata(FO.fig, FO); %reset the event counter
    
    set(correcttext, 'String',['Correct: ',num2str(length(correct))]);
    set(misstext, 'String', ['Miss: ',num2str(length(miss))]);
    set(FAtext, 'String', ['FA: ',num2str(length(falsealarm))]);
    set(countetext, 'String', ['Windows Remaining: ',num2str(counter)]);
end
    
%What to do after the loop?
%%
%Calculate number of misses/s, false alarms/s, corrects/s?
%Estimate miss% from miss./(correct+miss);
%Estimate FA % from FA./(correct+FA)
nummiss = length(miss);
numcorrect = length(correct);
numFA = length(falsealarm);
estmissperc = nummiss./(numcorrect+nummiss);
estFAperc = numFA./(numcorrect+numFA);

%After looking at Xmin of Y total detection minutes, we've found that your
%detector is missing an estimated 12% of events and falsely detecting 0%

EEoutput.ROCresults.FA = estFAperc;
EEoutput.ROCresults.TP = 1-estmissperc;
%%
close(FO.fig)
end

