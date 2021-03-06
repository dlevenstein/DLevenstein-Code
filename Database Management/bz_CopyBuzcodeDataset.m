function [  ] = bz_CopyBuzcodeDataset( fromDir,toDir )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%

% -a, --archive               archive mode; equals -rlptgoD (no -H,-A,-X)
% -u, --update                skip files that are newer on the receiver
% -v, --verbose               increase verbosity
% -z, --compress              compress file data during the transfer
% -L, --copy-links            transform symlink into referent file/dir 
% -P                          show progress during transfer
% -m                          no empty folders



%Include only buzcode things
unix(['rsync -mauvzLP ',...
    fromDir,' ',toDir,...
    ' --include=*.lfp --include=*.cellinfo.mat --include=*.states.mat',...
    ' --include=*.LFP.mat --include=*.events.mat --include=*.sessionInfo.mat',...
    ' --include=*.behavior.mat --include=*.popinfo.mat --include=*.eeg',...
    ' --include=*.xml --include=*.jpg --include=*.pdf',...
    ' --include=*/ --exclude=*'])  %exclude everything except:])

%Exclude .dat, .clu .res
% unix(['rsync -mauvzLP ',...
%     fromDir,' ',toDir,...
%     ' --exclude=*.dat --exclude=*.fet.* --exclude=*.res.*',...
%     ' --exclude=*.clu.* --exclude=*.spk.* --exclude=*.npy',...
%     ' --exclude=.phy'])

end

