function [ output_args ] = RunAnalysis(analysisfunction,datasetfolder)
%RunAnalysis(analysisfunction,datasetfolder) is a high level function that
%takes as its input an analysis function and returns the result of the
%analysis on multiple recordings
%
%INPUT
%   analysisfunction    A function name.  The function should be structured as
%                       [resultargs] = analysisfunction(datasetfolder,recordingname,figfolder)
%                       where resultargs is a variable number of analysis
%                       results.
%   datasetfolder       The folder in which the dataset lives. This
%                       requirement will be taken out in future versions.
%
%
%TO ADD
%   option to select from spreadsheet of recordings with recording
%   properties
%
%DLevenstein 2016
%% Select 
foldercontents = dir(datasetfolder);
possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
[s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                'ListString',possiblerecordingnames);
recordingname = possiblerecordingnames(s);

%% Make a figure folder
functionpath = which(analysisfunction);
functionpath = fileparts(functionpath);
figfolder = [fullfile(functionpath,'AnalysisFigs',analysisfunction),'/'];
if ~exist(figfolder,'dir')
    mkdir(figfolder)
end

%% Loop through the recordings and run the analysis
numresults = nargout(analysisfunction);
numrecs = length(recordingname);
display(['Running Analysis on Recordings (',num2str(numrecs),')'])
for rr = 1:numrecs
    [results{1:numresults}] = feval(analysisfunction,datasetfolder,recordingname{rr},figfolder);
    close all
end




end

