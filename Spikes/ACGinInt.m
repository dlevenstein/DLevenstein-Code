function [ACG,lags] = ACGinInt( spiketimes,int,binsize )
%UNTITLED Summary of this function goes here
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
dt = binsize;
[spikemat,t] = SpktToSpkmat(spiketimes,[],dt);
[~,int_ts,~] = RestrictInts(t,int);
spikemat(~int_ts,:) = nan;
%%
chopfront = find(~isnan(spikemat(:,1)),1,'first');
chopend = find(~isnan(spikemat(:,1)),1,'last');
spikemat = spikemat(chopfront:chopend,:);
spikemat(isnan(spikemat))=0;
%%
support = 0.5;
ACG = zeros((2*support)/dt+1,numcells);

for cc = 1:numcells
    [ACG(:,cc),lags] = xcorr(spikemat(:,cc),support/dt);
end

%%
lags = lags*binsize;
ACG(end/2+0.5,:) = 0;
%%
% figure
%     for cc = 1:numcells
%         bar(lags,ACG(:,cc))
%         pause
%     end
end







%Alternative to using xcorr: limit spiketimes added to isi to those within
%int/support, should make this go a lot faster.
%     function autocorrelogram(spiketimes)
%     %spiketimes = [10413177585,10413282812,10413379677,10413402313,10413410739,10413422026,..]
%         isi=[];
%         for i = 1:length(spiketimes)
%             isi = [isi; spiketimes-spiketimes(i)];
%         end
%         isi(isi>2000000)=[];
%         isi(isi<-2000000)=[];
%         isi=isi./1000000;
%         save isi;
%         figure
%         hist(isi,-2:0.02:2);

