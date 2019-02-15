function [ results,baseNames ] = GetMatResults(lookatfolder,matname,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   lookatfolder    folder in which the results of analysis are
%   matname         name of the analysis mat (i.e. '_NREMStats')
%
%   (options)
%       recordings      'select' to choose whichrecordings to load
%
%DLevenstein 2016
%% DEV
% lookatfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisScripts/AnalysisFigs/NREMTemporalStatsAnalysis';
% matname = '_NREMStats';
%%
p = inputParser;
addParameter(p,'select',false,@islogical)
parse(p,varargin{:})
select = p.Results.select;
%%
matnames = dir(fullfile(lookatfolder,['*',matname,'.mat']));
baseNames = cellfun(@(X) strtok(X,'.'),{matnames.name},'UniformOutput',false);

%User select recordings
if select
    
      [s,v] = listdlg('PromptString','Which recording(s) would you like?',...
                    'ListString',baseNames);
    matnames = matnames(s);
    baseNames = baseNames(s);
end

nummats = length(matnames);

FIELDMISMATCH=false;
%%
for mm = 1:nummats
    filename = fullfile(lookatfolder,matnames(mm).name); 
    loadmat = load(filename);
    
    %
    [pathstr,name] = fileparts(filename);
    loadmat.name = string(name);
    
    %Check if the new .mat has any additional fields
    if exist('results','var')
            [results,loadmat] = bz_Matchfields(results,loadmat,'remove');
    end  
    
%     %Check if the new .mat has any additional fields
%     matfields = fieldnames(loadmat);
%     resultsfields = fieldnames(results);
%     newfields = setdiff(matfields,resultsfields);
%     if ~isempty(newfields);
%         for ff = 1:length(newfields)
%             results(1).(newfields{ff}) = []; 
%             %results(1).name = name;
%         end
%         FIELDMISMATCH=true;
%     end
%     
%     %results(mm).name = name;
%     results = orderfields(results,loadmat);
    results(mm) = loadmat;
    %results(mm).name = name;
       
end

if FIELDMISMATCH
    display('One or more of your .mats has missing fields')
end

end

