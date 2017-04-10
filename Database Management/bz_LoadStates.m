function [ states ] = bz_LoadStates( statesName,baseName,datasetPath )
%[ states ] = bz_LoadStates( statesName,baseName,datasetPath ) function for
%loading states.mat files. states.mat files are saved as...
% datasetPath/baseName/baseName.statesName.states.mat
%
%statesName can be the name of a states.mat file, or can be 'all' to load
%all states.mat files for a given recording.
%
%DLevenstein 2017
%%
% if strcmp('statesName','all')
%     allStates = dir(fullfile(datasetPath,baseName,[baseName,'.','*','.states.mat']));
%     for


statesfile = fullfile(datasetPath,baseName,[baseName,'.',statesName,'.states.mat']);
evalin('caller',['load(''',statesfile,''')']);


end

