function [ paircorr ] = PairCorrinInt(spiketimes,int,binsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%%

%% Deal with Input Types
numcells = length(spiketimes);
if numcells == 0
    display('case no cells needs outputs')
    pause
    return
end

%Spiketimes can be: tsdArray of cells, cell array of cells, cell array of
%tsdArrays (multiple populations)
if isa(spiketimes,'tsdArray')
    numcells = length(spiketimes);
    for cc = 1:numcells
        spiketimestemp{cc} = Range(spiketimes{cc},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
elseif isa(spiketimes,'cell') && isa(spiketimes{1},'tsdArray')
    numpop = length(spiketimes);
    lastpopnum = 0;
    for pp = 1:numpop
        if length(spiketimes{pp})==0
            spiketimes{pp} = {};
            popcellind{pp} = [];
            continue
        end
        for cc = 1:length(spiketimes{pp})
            spiketimestemp{cc} = Range(spiketimes{pp}{cc},'s');
        end
        spiketimes{pp} = spiketimestemp;
        popcellind{pp} = [1:length(spiketimes{pp})]+lastpopnum;
        lastpopnum = popcellind{pp}(end);
        clear spiketimestemp
    end
    spiketimes = cat(2,spiketimes{:});
    numcells = length(spiketimes);
    subpop = 'done';
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


%%
dt = binsize./4;
overlap = 4;
[spikemat,t] = SpktToSpkmat(spiketimes,[],dt,overlap);
[~,int_ts,~] = RestrictInts(t,int);
spikemat(~int_ts,:) = [];

%% 

paircorr = corr(spikemat,'type','spearman');
paircorr(diag(diag(true(size(paircorr))))) = 0;
%%
% figure
%     imagesc(paircorr)
end

