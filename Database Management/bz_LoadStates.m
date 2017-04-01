function [ states ] = bz_LoadStates( statesName,baseName,datasetPath )
%[ states ] = bz_LoadStates( statesName,baseName,datasetPath ) function for
%loading states.mat files.
%
%
%DLevenstein 2017
%%

statesfile = fullfile(datasetPath,baseName,[baseName,'.',statesName,'.states.mat']);

evalin('base',['load(''',statesfile,''')']);

%assignin('base',statesName

end

