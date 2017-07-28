function [ output_args ] = EventExplorer(eventsname,basePath )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
eventsname = 'SlowWave';

%%
if ~exist('basePath','var')
    basePath = pwd;
end
baseName = bz_BasenameFromBasepath(basePath);

%%
%eventsfile = fullfile(basePath,[baseName,'.',eventsname,'.events.mat']);

%events = bz_LoadEvents(eventsname);

events = bz_LoadStates(basePath,eventsname);
eventstype = 'states';

switch eventstype
    case 'states'
        intnames = fieldnames(events.ints);
        [s,v] = listdlg('PromptString','Which interval type to explore?',...
            'SelectionMode','single','ListString',intnames);
        exploreint = intnames{s};
end
       
end

