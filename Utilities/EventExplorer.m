function [ output_args ] = EventExplorer(eventsname,basePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%DLevenstein 2017
%%
eventsname = 'SlowWave';
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

%events = bz_LoadEvents(eventsname);

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
FO.fig = figure;
set(FO.fig, 'numbertitle', 'off', 'name', ['Recording: ', FO.baseName,'. Events: ',FO.EventName]);
set(FO.fig, 'Tag', 'EventExplorerMaster');

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
numevents = 50; %number of events to look at. determine to maximize sampling or have user input (FO.uinput.field?)
randevents = randi(length(FO.EventTimes),[numevents,1]);
randevents = FO.EventTimes(randevents);

%Find any events that are within winsize of another event
closeevents = abs(diff(randevents))<winsize;
% numtries =0;
% while any(closeevents) && numtries<50
%     randevents(closeevents)
% end
%Problem:overlapping ints
%[ ints ] = MergeSeparatedInts( ints,minseparation )

%% Loop figure to show events and get user input

counter = numevents;
miss=[];correct=[];falsealarm=[];
QUITLOOP = false;
while counter>0 && QUITLOOP~=true
    %EventUI(FO)   %function that will run the user interface
    thiseventtime = randevents(counter);
    thiseventwin = thiseventtime+winsize.*[-0.5 0.5];
    inwinevents = FO.EventTimes(FO.EventTimes>=thiseventwin(1) &FO.EventTimes<=thiseventwin(2));
    
    lfpwin = subplot(3,1,2);
    bz_MultiLFPPlot(lfp,'timewin',thiseventwin,'spikes',spikes)
    hold on
    %plot(spikes.spindices(:,1),spikes.spindices(:,2),'.')
    plot(inwinevents,zeros(size(inwinevents)),'r+')
    
    %Get user feedback
    title({'Left click missed events. Right click incorrectly ID''d events.',...
        'RETURN when complete. ESC,RETURN to finish early.',...
        [num2str(counter),' Windows Remain']})
    [x,y,button] = ginput;
    miss = [miss; x(button==1)];
    falsealarm = [falsealarm; interp1(inwinevents,inwinevents,x(button==3),'nearest')]; 
    correct = [correct; inwinevents(~ismember(inwinevents,falsealarm))];
    if any(button==27);QUITLOOP=true;end
    
    %Update window to show how many misses/FAs have been ID'd
    %When click.... mark the event as FA or miss on screen...
    %Top window to give context?
    %Instructions window
    
    hold off
    counter = counter-1;
    
end
    
%%
%Calculate number of misses/s, false alarms/s, corrects/s?
%Estimate miss% from miss./(correct+miss);
%Estimate FA % from FA./(correct+FA)
nummiss = length(miss);
numcorrect = length(correct);
numFA = length(falsealarm);
estmissperc = nummiss./(numcorrect+nummiss)
estFAperc = numFA./(numcorrect+numFA)

%After looking at Xmin of Y total detection minutes, we've found that your
%detector is missing an estimated 12% of events and falsely detecting 0%

end

