function [ behavior ] = bz_LoadBehavior( behaviorName,baseName,datasetPath )
%[ states ] = bz_LoadStates( statesName,baseName,datasetPath ) function for
%loading states.mat files.
%
%
%DLevenstein 2017
%%

behaviorfile = fullfile(datasetPath,baseName,[baseName,'.',behaviorName,'.behavior.mat']);

evalin('base',['load(''',behaviorfile,''')']);

%assignin('base',statesName

end

