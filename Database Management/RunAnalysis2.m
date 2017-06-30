function [  ] = RunAnalysis(analysisfunction,lookatfolder)
%RunAnalysis(analysisfunction,datasetfolder) is a high level function that
%takes as its input an analysis function and returns the result of the
%analysis on multiple recordings
%
%INPUT
%   analysisfunction    A function name.  The function should be structured as
%                       [resultargs] = analysisfunction(basePath,figfolder)
%                       where resultargs is a variable number of analysis
%                       results.
%   lookatfolder       (optional) The folder in which the dataset lives. This
%                       requirement will be taken out in future versions.
%                       If not included, will consult datasetguide.xlsx
%
%
%TO ADD
%   option to select from spreadsheet of recordings with recording
%   properties
%
%DLevenstein 2017
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
%         dropboxdesktop = '/mnt/data1/Dropbox/research/';
%         dropboxlaptop = '/Users/dlevenstein/Dropbox/Research/';
%         if exist(dropboxdesktop)
%             dropboxdatabasefolder = fullfile(dropboxdesktop,'/Datasets');
%         elseif exist(dropboxlaptop)
%             dropboxdatabasefolder = fullfile(dropboxlaptop,'/Datasets');
%         else display('No dropbox folder!'); end
%    % dropboxdatabasefolder = '/Users/dlevenstein/Dropbox/Research/Datasets';

    dropboxdatabasefolder = '/mnt/proraidDL/ProjectDatasets/S1State/'; %fix later
    datasetguidefilename = fullfile(dropboxdatabasefolder,'DatasetGuide.xlsx');
    
    datasetguide=readtable(datasetguidefilename);
    possiblebaseNames = datasetguide.RecordingName;
    possiblebasePaths = datasetguide.basePath;
    possibleGenotypes = datasetguide.Genotype; %make this general later
 %   possibledesktopfolders = datasetguide.DesktopFolder;
    [s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                    'ListString',possiblebaseNames);
    baseName = possiblebaseNames(s);
%     basePath = cellfun(@(X) fullfile(dropboxdatabasefolder,X),...
%         possiblebasePaths(s),'UniformOutput',false);
    basePath = possiblebasePaths(s); 
    genotypes = possibleGenotypes(s);

    %%Select Recordings from Dataset Folder 
    case 'datasetfolder'
%     foldercontents = dir(datasetfolder);
%     possiblerecordingnames = {foldercontents([foldercontents.isdir]==1).name};
    [possiblebasePaths,possiblebaseNames] = bz_FindBasePaths(datasetfolder);
    [s,v] = listdlg('PromptString','Which recording(s) would you like analyze?',...
                    'ListString',possiblebaseNames);
    baseName = possiblebaseNames(s);
    basePath = possiblebasePaths(s); 
    [temp{1:length(s)}] = deal(datasetfolder);
    datasetfolder = temp; 
    clear temp
    
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
numrecs = length(baseName);
display(['Running Analysis on Recordings (',num2str(numrecs),')'])
for rr = 1:numrecs
    display(['Recording: ',num2str(rr),' of ',num2str(numrecs),' - ',baseName{rr}])
    
    [inputNames, outputNames] = get_arg_names(fullfile(functionpath,[analysisfunction,'.m']));
    outputNames = outputNames{1}; inputNames = inputNames{1};
    

        %[results{1:numresults}] = feval(analysisfunction,datasetfolder{rr},recordingname{rr},figfolder);
        [results{1:numresults}] = feval(analysisfunction,basePath{rr},figfolder);
        
    for ss = 1:numresults
        eval([outputNames{ss} '= results{ss};']);
    end
    
    recname = baseName{rr};
    %genotype = genotypes{rr};
    save([figfolder,baseName{rr},'_',analysisfunction],...
    outputNames{:});%,'genotype','recname')
    
    close all
end


end

