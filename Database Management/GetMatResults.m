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

FIELDMISMATCH=false;
%%
for mm = 1:nummats
    matname = fullfile(lookatfolder,matnames(mm).name); 
    loadmat = load(matname);
    
    %For the first .mat
    if ~exist('results','var'); results(mm) = loadmat; continue; end  
    
    %Check if the new .mat has any additional fields
    matfields = fieldnames(loadmat);
    resultsfields = fieldnames(results);
    newfields = setdiff(matfields,resultsfields);
    if ~isempty(newfields);
        for ff = 1:length(newfields)
            results(1).(newfields{ff}) = []; 
        end
        FIELDMISMATCH=true;
    end
    
    results = orderfields(results,loadmat);
    results(mm) = loadmat;
       
end

if FIELDMISMATCH
    display('One or more of your .mats has missing fields')
end

end

