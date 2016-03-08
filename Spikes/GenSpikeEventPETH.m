function [ ratePETH,popsynchPETH,popcellind ] = GenSpikeEventPETH(spiketimes,int,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   spiketimes
%   int
%   (optional)
%       'normwin'   %Sorrounding time to add to norm (in normalized time)
%       'synchdt'   %Synchrony dt                         (default: 0.005s)
%       'synchwin'  %Synchrony time window                 (default: 0.04s)
%       'twin'      %Time window around onset/offset PETH (default: [2 2]s)
%       'dt_onoff'  %binsize for onset/offset PETH         (default: 0.05s)
%       'normbins'  %Number of bins in the time normalized interval 0-1
%       'sorttype'
%       'subpop'
%
%OUTPUT
%   ratePETH
%       .onset
%       .offset
%       .norm
%       .t_on
%       .t_off
%       .t_norm
%       .cellpopidx
%   popsynchPETH
%       .onset
%       .offset
%       .norm
%       .t_on
%       .t_off
%       .t_norm
%
%TO DO
%   -...lots

%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultNormwin = 1;    %Time overhang before/after to add (in normalized time)

defaultSynchdt = 0.005;
defaultSynchwin = 0.04;

defaultTwin = [2 2];   %[Before/after onset/offset, after/before onset/offset]
defaultDt_onoff = 0.05; %dt bin for onset/offset aligned

defaultNormbins = 30;  %Number of bins in the time normalized interval 0-1

defaultSorttype = 'rate';
validSorttypes = {'pca','none','sortf','sortorder','celltype','rate'};
checkSorttype = @(x) any(validatestring(x,validFranges)) || size(x) == [1,2];

defaultDOWN = false;

addParameter(p,'normwin',defaultNormwin,@isnumeric)
addParameter(p,'synchdt',defaultSynchdt,@isnumeric)
addParameter(p,'synchwin',defaultSynchwin,@isnumeric)
addParameter(p,'twin',defaultTwin,@isnumeric)
addParameter(p,'dt_onoff',defaultDt_onoff,@isnumeric)
addParameter(p,'normbins',defaultNormbins,@isnumeric)
addParameter(p,'sorttype',defaultSorttype,checkSorttype)
addParameter(p,'subpop',0,@isnumeric)

parse(p,varargin{:})
%Clean up this junk...
normwin = p.Results.normwin; 
synchdt = p.Results.synchdt;
synchwin = p.Results.synchwin;
subpop = p.Results.subpop;
twin = p.Results.twin;
dt_onoff = p.Results.dt_onoff;
normbins = p.Results.normbins;


%% Time parms
totnormtime = 2*normwin+1;
totnormbins = totnormtime*normbins;


%% Deal with input types

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


%% Subpopulations
switch subpop
    case 0
        numpop = 1;
        popcellind = {1:numcells};
    case 'done'
    otherwise
        pops = unique(subpop);
        numpop = length(pops);
        for pp = 1:numpop
            popcellind{pp} = find(subpop == pops(pp));
        end
end

%% Ints for time norm and onset/offsets

intlen = int(:,2)-int(:,1);
windur = twin(2)+twin(1);

%Norm ints
normints = [int(:,1)-normwin*intlen int(:,2)+normwin*intlen];
%ratemat.t_norm = linspace(-normwin,normwin+1,totnormbins);

%Onset ints
onints = [int(:,1)-twin(1) int(:,1)+twin(2)];
%t_on = [-twin(1):dt_onoff:twin(2)];

%Offset ints
offints = [int(:,2)-twin(2) int(:,2)+twin(1)];
%t_off = [-twin(2):dt_onoff:twin(1)];

%% Time Norm Spikes

for cc = 1:numcells
    [timenormspiketimes{cc},~,~,overlaps] = SortByIntTime(spiketimes{cc},normints,'norm');
    timenormspiketimes{cc} = [timenormspiketimes{cc};[overlaps{:}]'];
end
dt = 1/totnormbins;
[ratePETH.norm,ratePETH.t_norm] = SpktToSpkmat(timenormspiketimes,[0 1],dt);
ratePETH.t_norm = ratePETH.t_norm*totnormtime-normwin;
%% Onset/Offset Spikes
for cc = 1:numcells
    [onsetspiketimes{cc},~,~,overlaps] = SortByIntTime(spiketimes{cc},onints,'onset');
    onsetspiketimes{cc} = [onsetspiketimes{cc};[overlaps{:}]'];
end
[onsetspikemat,t_on] = SpktToSpkmat(onsetspiketimes,[0 windur],dt_onoff);
ratePETH.onset = (onsetspikemat./length(onints(:,1)))./dt_onoff;
ratePETH.t_on = t_on-twin(1);

for cc = 1:numcells
    [offsetspiketimes{cc},~,~,overlaps] = SortByIntTime(spiketimes{cc},offints,'onset');
    offsetspiketimes{cc} = [offsetspiketimes{cc};[overlaps{:}]'];
end
[offsetspikemat,t_off] = SpktToSpkmat(offsetspiketimes,[0 windur],dt_onoff);
ratePETH.offset = (offsetspikemat./length(offints(:,1)))./dt_onoff;
ratePETH.t_off = t_off-twin(2);

ratePETH.cellpopidx = zeros(1,numcells);
%% Pop Synch

% Calculate Population Synchrony Coupling
overlap = synchwin/synchdt;
[spikemat,t_synch] = SpktToSpkmat(spiketimes, [], synchdt,overlap);


for pp = 1:numpop
    if length(popcellind{pp}) == 0

        popsynchPETH(pp).t_norm = ratePETH.t_norm;
        popsynchPETH(pp).norm = nan(size(popsynchPETH(pp).t_norm));
        popsynchPETH(pp).t_on = [-twin(1):synchdt:twin(2)]';
        popsynchPETH(pp).onset = nan(size(popsynchPETH(pp).t_on));
        popsynchPETH(pp).t_off =[-twin(2):synchdt:twin(1)]';
        popsynchPETH(pp).offset = nan(size(popsynchPETH(pp).t_off));
        numpop = numpop-1;
        continue
    end
    %%
    ratePETH.cellpopidx(popcellind{pp}) = pp;
    numpopcells = length(popcellind{pp});
    popsynch = sum(spikemat(:,popcellind{pp})>0,2)./numpopcells;

    %% Time Norm
    popsynch_normepochs = IsolateEpochs2(popsynch,normints,0,1/synchdt);
    popsynch_normepochs = TimeNormalize(popsynch_normepochs,totnormbins);
    popsynchPETH(pp).norm = mean(cat(2,popsynch_normepochs{:}),2);
    popsynchPETH(pp).t_norm = ratePETH.t_norm;
    %% Onset
    popsynch_onepochs = IsolateEpochs2(popsynch,onints,0,1/synchdt);
    popsynchPETH(pp).onset = mean(cat(2,popsynch_onepochs{:}),2);
    popsynchPETH(pp).t_on = [-twin(1):synchdt:twin(2)]';
    
    %% Onset
    popsynch_offepochs = IsolateEpochs2(popsynch,offints,0,1/synchdt);
    popsynchPETH(pp).offset = mean(cat(2,popsynch_offepochs{:}),2);
    popsynchPETH(pp).t_off = [-twin(2):synchdt:twin(1)]';
end

%% Figure
figure
    subplot(6,1,1:2)
        imagesc(ratePETH.t_norm,1:numcells,(ratePETH.norm)')
        %colorbar
        xlim([-0.8 1.8])
    subplot(6,1,3)
        hold on
        for pp = 1:numpop
            plot(ratePETH.t_norm,popsynchPETH(pp).norm)
        end
        xlim([-0.8 1.8])
    subplot(6,2,[7,9])
        imagesc(ratePETH.t_on,1:numcells,(ratePETH.onset)')
        %colorbar
        xlim([-twin(1) twin(2)])
    subplot(6,2,11)
        hold on
        for pp = 1:numpop
            plot(popsynchPETH(pp).t_on,popsynchPETH(pp).onset)
        end
        xlim([-twin(1) twin(2)])
    subplot(6,2,[8,10])
        imagesc(ratePETH.t_off,1:numcells,(ratePETH.offset)')
        %colorbar
        xlim([-twin(2) twin(1)])
    subplot(6,2,12)
        hold on
        for pp = 1:numpop
            plot(popsynchPETH(pp).t_off,popsynchPETH(pp).offset)
        end
        xlim([-twin(2) twin(1)])
