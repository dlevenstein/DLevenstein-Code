function [ EEoutput ] = EventExplorer(events,basePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   events      Name of events in a string (i.e. 'SlowWaves') or buzcode 
%               events structure
%   basePath    default: pwd
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

if isstring(events) || ischar(events)
    eventsname = events;
    [events,FO.eventsfilename] = bz_LoadStates(basePath,events);
    eventstype = 'states';

    switch eventstype
        case 'states'
            intnames = [fieldnames(events.ints);fieldnames(events)]; %remove second term here
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
FO.EventName = eventsname;
FO.basePath = basePath;
%% Make the EventExplorer Figure
%EE_Initiate( FO )
usechans = events.detectorinfo.detectionparms.SWchannel;
[ SleepState ] = bz_LoadStates(FO.basePath,'SleepState');
restrictint = SleepState.ints.NREMstate;
%lfp = bz_GetLFP(usechans,'basepath',FO.basePath,'intervals',restrictint);
FO.data.lfp = bz_GetLFP(usechans,'basepath',FO.basePath);
FO.data.spikes = bz_GetSpikes('basepath',FO.basePath);
%%

%Position for the main interface
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
set(FO.fig, 'menubar', 'none');
%set(FO.fig, 'CloseRequestFcn', {@CloseDialog}); The close function
set(FO.fig,'WindowButtonDownFcn', {@MouseClick});

%From StateEditor
% FO.fig = figure('KeyReleaseFcn', {@DefKey}, 'Position', posvar);
% set(FO.fig, 'numbertitle', 'off', 'name', ['Event: ', FO.baseName]);
% set(FO.fig,'WindowButtonDownFcn', {@MouseClick}, 'WindowButtonUpFcn', {@unMouseClick}, 'Units', 'normalized');
% set(FO.fig, 'WindowButtonMotionFcn', {@Nothing}, 'WindowScrollWheelFcn', {@MouseScroll});
% set(FO.fig, 'CloseRequestFcn', {@CloseDialog});
% set(FO.fig, 'Tag', 'EventExplorerMaster');

%Set up the view window
FO.viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%Event Selection panel

FO.winsize = 8;
FO.currevent = 1;
EventVewPlot(FO)

%Set up the navigation panel
FO.NavPanel = uipanel('FontSize',12,...
        'Position',[.65 .20 0.25 0.15]);
nextbtn = uicontrol('Parent',FO.NavPanel,...
    'Position',[160 70 100 40],'String','->',...
     'Callback','FO.currevent=FO.currevent+1;EventVewPlot(FO)');
prevbtn = uicontrol('Parent',FO.NavPanel,...
    'Position',[50 70 100 40],'String','<-',...
     'Callback','FO.currevent=FO.currevent-1;EventVewPlot(FO)');
 editwinsize  = uicontrol('Parent',FO.NavPanel,'Style','edit',...
    'Position',[130 20 60 25],'String',num2str(FO.winsize),...
    'Callback','FO.winsize=str2num(editwinsize.String);EventVewPlot(FO)');
winsizetext = uicontrol('Parent',FO.NavPanel,...
    'Position',[50 10 80 30],'style','text',...
    'string','Window Size (s):','HorizontalAlignment','left'); 

%Process selection panel (i.e. event detection, other?)
    
%Store the data in the figure - do this at the end of each function?
guidata(FO.fig, FO);


%% THIS IS FOR THE EVENT DETECIONT REVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function: DetectionReview(FO)
%% Select the time windows to look at 
%User input for this
numevents = 5; %number of events to look at. determine to maximize sampling or have user input (FO.uinput.field?)
%Selecting from events:
%randevents = randsample(FO.EventTimes,numevents);
%Selecting from random times
[NREMtime,~,~] = RestrictInts(FO.data.lfp.timestamps,restrictint);
randevents = randsample(NREMtime,numevents);
%Find any events that are within winsize of another event to remove them?
%Also should look at windows in which no events were detected?  Maybe just
%random windows in the detectionwin
closeevents = abs(diff(randevents))<FO.winsize;
% numtries =0;
% while any(closeevents) && numtries<50
%     randevents(closeevents)
% end
%Problem:overlapping ints
%[ ints ] = MergeSeparatedInts( ints,minseparation )


%%
%Store the current user action
FO.currentuseraction = 'MarkEvents';
FO.markedevents = []; %clear any previously marked events
guidata(FO.fig, FO); %store the data in the figure

counter = numevents;
miss=[];hit=[];falsealarm=[];
lookedatwins = [];

%UI panel for event review
FO.EventPanel = uipanel('Title','Detection Review','FontSize',12,...
        'Position',[.65 .05 0.25 0.3]);
%Buttons for next loop or finishing
nextbtn = uicontrol('Parent',FO.EventPanel,...
    'Position',[160 20 125 40],'String','Next Window (Return)',...
          'Callback','uiresume(gcbf)');
qutbtn = uicontrol('Parent',FO.EventPanel,...
    'Position',[160 65 125 40],'String','Quit Early (Esc)',...
         'Callback','QUITLOOP=true;uiresume(gcbf)');
%Instruction text
instructtext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 100 400 100],'style','text',...
    'string',{'Instructions:',...
    '   - LEFT click missed events (o)','   - RIGHT click False Alarms (x)'},...
    'HorizontalAlignment','left'); 
%Display counter text
correcttext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 60 125 15],'style','text',...
    'string','Correct: 0','HorizontalAlignment','left'); 
misstext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 40 125 15],'style','text',...
    'string','Miss: 0','HorizontalAlignment','left');
FAtext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 20 125 15],'style','text',...
    'string','FA: 0','HorizontalAlignment','left');
countetext = uicontrol('Parent',FO.EventPanel,...
    'Position',[25 100 125 15],'style','text',...
    'string',['Windows Remaining: ',num2str(counter)],'HorizontalAlignment','left');

%The Event Review Loop
QUITLOOP = false;
while counter>0 && QUITLOOP~=true
    %EventUI(FO)   %function that will run the user interface
    thiseventtime = randevents(counter);  
    viewinfo = EventVewPlot( FO,thiseventtime);
    uiwait(FO.fig)  %Wait here until the user clicks quit or next
    
    lookedatwins = [lookedatwins; viewinfo.thiseventwin];
    inwinevents = viewinfo.inwinevents;
    
    
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
    hit = [hit; inwinevents(~ismember(inwinevents,falsealarm))];
    FO.markedevents = [];  guidata(FO.fig, FO); %reset the marked events
    
    counter = counter-1;
    set(correcttext, 'String',['Correct: ',num2str(length(hit))]);
    set(misstext, 'String', ['Miss: ',num2str(length(miss))]);
    set(FAtext, 'String', ['FA: ',num2str(length(falsealarm))]);
    set(countetext, 'String', ['Windows Remaining: ',num2str(counter)]);
end
%Calculate total number of miss,hit,FA
numMiss = length(miss);
numHit = length(hit);
numFA = length(falsealarm);    

%Calculate total amount of time/percentage of detection time (detectionintervals) looked at

%Put things in the output structure
EEoutput.EventReview.lookedatwins = lookedatwins;
EEoutput.EventReview.miss = miss; 
EEoutput.EventReview.hit = hit;
EEoutput.EventReview.falsealarm = falsealarm;
EEoutput.EventReview.estMissperc = numMiss./(numHit+numMiss);
EEoutput.EventReview.estFAperc = numFA./(numHit+numFA);
EEoutput.EventReview.ReviewDate = today;
EEoutput.EventReview.EventsType = FO.EventName;

%UI: Done!  Would you like to save the results to (eventsfilename?)
%Make function that does this: SaveResults(FO,EEoutput)
if isfield(FO,'eventsfilename')
    button = questdlg(['Event Review Complete! Would you like to add the results to ',...
        FO.eventsfilename,'?'],'Good Job!');
    switch button
        case 'Yes'
            %Load the events file, add the field, save the events file
            try %Only do this if the correct named structure lives in the file
                eventsfile = load(FO.eventsfilename,FO.EventName);
                eventsfile.(FO.EventName).EventReview = EEoutput.EventReview;
                save(FO.eventsfilename,'-struct','eventsfile',FO.EventName,'-append')
            catch
                display([' Save failed... ',FO.eventsfilename,' may not ',...
                    'contain a structure titled ',FO.EventName,'.',])
            end
    end
end
%%
close(FO.fig)
end

