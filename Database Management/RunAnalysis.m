function [ output_args ] = RunAnalysis(analysisfunction,lookatfolder)
%RunAnalysis(analysisfunction,datasetfolder) is a high level function that
%takes as its input an analysis function and returns the result of the
%analysis on multiple recordings
%
%INPUT
%   analysisfunction    A function name.  The function should be structured as
%                       [resultargs] = analysisfunction(datasetfolder,recordingname,figfolder)
%                       where resultargs is a variable number of analysis
%                       results.
%                       (optional) function can be structured as [resultargs] =
%                        analysisfunction(datasetfolder,recordingname,figfolder,desktopfolder)
%                       if desktopfolder is needed
%   lookatfolder       (optional) The folder in which the dataset lives. This
%                       requirement will be taken out in future versions.
%                       If not included, will consult datasetguide.xlsx
%
%
%TO ADD
%   option to select from spreadsheet of recordings with recording
%   properties
%
%DLevenstein 2016
%% Select Recordings to Analyze 
if ~exist('lookatfolder','var')
    selectionmode = 'datasetguide';
else
    selectionmode = 'datasetfolder';
    datasetfolder = lookatfolder;
end

switch selectionmode
              
    %%Select Recordings from DatasetGuide spreadsheet
    case 'datasetguide'
        dropboxdesktop = '/mnt/data1/Dropbox/research/';
        dropboxlaptop = '/Users/dlevenstein/Dropbox/Research/';
        if exist(dropboxdesktop)
            dropboxdatabasefolder = fullfile(dropboxdesktop,'/Datasets');
        elseif exist(dropboxlaptop)
            dropboxdatabasefolder = fullfile(dropboxlaptop,'/Datasets');
        else display('No dropbox folder!'); end
   % dropboxdatabasefolder = '/Users/dlevenstein/Dropbox/Research/Datasets';
    datasetguidefilename = fullfile(dropboxdatabasefolder,'DatasetGuide.xlsx');
    datasetguide=readtable(datasetguidefilename);
    possiblerecordingnames = datasetguide.RecordingName;
    possibledatasetfolders = datasetguide.Dataset;
    possibledesktopfolders = datasetguide.DesktopFolder;
    [s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                    'ListString',possiblerecordingnames);
    recordingname = possiblerecordingnames(s);
    datasetfolder = cellfun(@(X) fullfile(dropboxdatabasefolder,X),...
        possibledatasetfolders(s),'UniformOutput',false);
    desktopfolder = possibledesktopfolders(s);
    

    %%Select Recordings from Dataset Folder 
    case 'datasetfolder'
    foldercontents = dir(datasetfolder);
    possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
    [s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                    'ListString',possiblerecordingnames);
    recordingname = possiblerecordingnames(s);
    
end

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
    display(['Recording: ',num2str(rr),' of ',num2str(numrecs),' - ',recordingname{rr}])
    
    [inputNames, outputNames] = get_arg_names(fullfile(functionpath,[analysisfunction,'.m']));
    outputNames = outputNames{1}; inputNames = inputNames{1};
    
    if length(inputNames)==4
        [results{1:numresults}] = feval(analysisfunction,datasetfolder{rr},recordingname{rr},figfolder,desktopfolder{rr});
    else
        [results{1:numresults}] = feval(analysisfunction,datasetfolder{rr},recordingname{rr},figfolder);
    end
        
    for ss = 1:numresults
        eval([outputNames{ss} '= results{ss};']);
    end
    
    %recname = recordingname{rr};
    save([figfolder,recordingname{rr},'_',analysisfunction],...
    outputNames{:})%,'recname')
    
    close all
end


end

