function [ results ] = GetMatResults(lookatfolder,matname,recordings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   lookatfolder    folder in which the results of analysis are
%   matname         name of the analysis mat (i.e. '_NREMStats')
%   recordings      (optional)
%
%DLevenstein 2016
%% DEV
% lookatfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisScripts/AnalysisFigs/NREMTemporalStatsAnalysis';
% matname = '_NREMStats';

%%
matnames = dir(fullfile(lookatfolder,['*',matname,'.mat']));
nummats = length(matnames);

%%
for mm = 1:nummats
    matname = fullfile(lookatfolder,matnames(mm).name);        
    results(mm) = load(matname);    
end

end

