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
events = 'SlowWaves';
basePath = '/Users/dlevenstein/Dropbox/Research/Datasets/20140526_277um';
%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);


%%

%A0 = EE_Initiate(
%% Select the Events
%eventsfile = fullfile(basePath,[baseName,'.',eventsname,'.events.mat']);

if isstring(events) || ischar(events)
    eventsname = events;
    [events,FO.eventsfilename] = bz_LoadEvents(basePath,events);
    eventstype = 'events';
    
    if isempty(events)
        display('no events.mat by that name, trying to find a states.mat')
        [events,FO.eventsfilename] = bz_LoadStates(basePath,events);
        eventstype = 'states';
    end
end

switch eventstype
    case 'events'
        exploreint = events.timestamps;
        exploreintname = events.detectorinfo.detectorname;

    case 'states'
        intnames = [fieldnames(events.ints)]; %remove second term here
        [s,v] = listdlg('PromptString','Which interval type to explore?',...
            'SelectionMode','single','ListString',intnames);
        exploreintname = intnames{s};
        try
        exploreint = events.ints.(exploreintname);
        catch
            exploreint = events.(exploreintname);
        end

end
%% 
FO.baseName = baseName;
FO.EventTimes = exploreint;
FO.EventName = eventsname;
FO.basePath = basePath;

%FAs and Misses if DetectionReview has already been run
REVIEWDONE = false;
if isfield(events,'EventReview')
    REVIEWDONE=true;
    FO.EventReview.miss = events.EventReview.miss;
    FO.EventReview.falsealarm = events.EventReview.falsealarm;
    FO.EventReview.miss(isnan(FO.EventReview.miss))=[];
    FO.EventReview.falsealarm(isnan(FO.EventReview.falsealarm))=[];
end
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
FO.fig = figure('KeyPressFcn', {@KeyDefinitions},'Position', posvar);
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

FO.currentuseraction = 'none';

%Set up the view window
FO.viewwin = subplot(3,1,2,'ButtonDownFcn',@MouseClick);
%Event Selection panel

FO.scaleLFP = 1;
FO.winsize = 8;
FO.currevent = 1;
FO.viewmode = 'event';


%Text of hotkey definitions

%Set up the navigation panel
FO.NavPanel = uipanel('FontSize',12,...
        'Position',[.65 .15 0.25 0.2]);
nextbtn = uicontrol('Parent',FO.NavPanel,...
    'Position',[160 70 100 40],'String','->',...
     'Callback',@NextEvent);
prevbtn = uicontrol('Parent',FO.NavPanel,...
    'Position',[50 70 100 40],'String','<-',...
     'Callback',@PrevEvent);
editwinsize  = uicontrol('Parent',FO.NavPanel,'Style','edit',...
    'Position',[170 20 60 25],'String',num2str(FO.winsize),...
    'Callback',@EditWinSize);
winsizetext = uicontrol('Parent',FO.NavPanel,...
    'Position',[50 10 120 30],'style','text',...
    'string','Window Size (s):','HorizontalAlignment','left'); 

if REVIEWDONE
    FO.eventtypeselection = uibuttongroup('Position',[0.1,0.1,0.1,0.15],'Visible','on',...
        'SelectionChangedFcn',@(bg,event) EventTypeSelector(bg,event));
    
    r1 = uicontrol(FO.eventtypeselection,'Style',...
                      'radiobutton',...
                      'String','event',...
                      'Position',[10 70 75 30]);

    r2 = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','misses',...
                      'Position',[10 40 75 30]);

    r3 = uicontrol(FO.eventtypeselection,'Style','radiobutton',...
                      'String','FAs',...
                      'Position',[10 10 75 30]);

end

    
%Process selection panel (i.e. event detection, other?)
    
%Store the data in the figure - do this at the end of each function?
guidata(FO.fig, FO);
EventVewPlot;
%uiwait(FO.fig) 

%[ EventReview ] = DetectionReview( )

end %Gen function end.


function NextEvent(hObject,eventdata)
    FO = guidata(hObject); 
    FO.currevent=FO.currevent+1;
    guidata(FO.fig, FO);
    EventVewPlot;
end

function PrevEvent(hObject,eventdata)
    FO = guidata(hObject); 
    FO.currevent=FO.currevent-1;
    guidata(FO.fig, FO);
    EventVewPlot;
end

function EditWinSize(hObject,eventdata)
    FO = guidata(hObject);
    FO.winsize=str2num(hObject.String);
    guidata(FO.fig, FO);
    EventVewPlot;
end

function EventTypeSelector(source,event)
    FO = guidata(source);
    FO.viewmode = event.NewValue.String;
    FO.currevent = 1;

    guidata(FO.fig, FO);
    EventVewPlot;
end