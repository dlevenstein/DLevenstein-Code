function [ states ] = bz_LoadStates( basePath,statesName )
%[ states ] = bz_LoadStates(basePath, statesName) function for
%loading states.mat files. states.mat files are saved as...
% datasetPath/baseName/baseName.statesName.states.mat
%
%statesName can be the name of a states.mat file, or can be 'all' to load
%all states.mat files for a given recording. If empty, prompts the user
%with a list of available states.mat files in basePath
%
%DLevenstein 2017
%%

[datasetPath,baseName] = fileparts(basePath);

% if strcmp('statesName','all')
%     allStates = dir(fullfile(datasetPath,baseName,[baseName,'.','*','.states.mat']));
%     for


statesfile = fullfile(datasetPath,baseName,[baseName,'.',statesName,'.states.mat']);
evalin('caller',['load(''',statesfile,''')']);


end

