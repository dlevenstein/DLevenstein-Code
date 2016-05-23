function [inttime4spike,phase4spike,...
    intphasehist,inthistbins,phasehistbins,...
    intlockbins,inttnumspks,intlockmag,intlockangle] = IntSpikePhase(spiketimes,LFP,frange,int,sf_LFP)
%IntSpikePhase( spiketimes,LFP,frange,int,sf_LFP ) gives the phase of
%spikes as a function of normalized time through a set of intervals (in
%progess)
%% DEV

% datasetfolder = '/Users/dlevenstein/Dropbox/Research/Datasets/BWData/';
% recname = '20140526_277um';
% %recname = 'Dino_061814_mPFC';
% figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisFigsMats/Spindle LFP/';
%
% load([datasetfolder,recname,'/',recname,'_LFP.mat'])
% load([datasetfolder,recname,'/',recname,'_StateIntervals.mat'])
% load([datasetfolder,recname,'/',recname,'_SSubtypes.mat'])
%
%
% spiketimes = Si;
% frange = [10 20];
% int = StateIntervals.Spindles;
% sf_LFP = 1250;
%% Optional Inputs

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

t_LFP = (1:length(LFP))'/sf_LFP;

%% Filter LFP
[~,~,LFP_phase] = FiltNPhase(LFP,frange,sf_LFP);

%% For each spike: phase and interval time
phase4spike = cellfun(@(X) interp1(t_LFP,LFP_phase,X,'nearest'),spiketimes,'UniformOutput',false);
inttime4spike = cellfun(@(X) SortByIntTime(X,int,'norm'),spiketimes,'UniformOutput',false);

%% Phase-interval spike histogram
numbins = 15;
phasehistbins = linspace(-pi,pi,numbins+1);
phasehistbins = phasehistbins(1:end-1)+0.5*diff(phasehistbins([1 2]));
inthistbins = linspace(0,1,numbins+1);
inthistbins = inthistbins(1:end-1)+0.5*diff(inthistbins([1 2]));
intphasehist = cellfun(@(X,Y) hist3([X,Y],{inthistbins,phasehistbins}),inttime4spike,phase4spike,'UniformOutput',false);
intphasehist = cat(3,intphasehist{:});

%% Phase Locking over int time

numlockbins = 5;
intlockbins = linspace(0,1,numlockbins+1);
[inttnumspks,~,binintt] = cellfun(@(X) histcounts(X,intlockbins),inttime4spike,'UniformOutput',false);
inttnumspks = cat(1,inttnumspks{:})';
intlockbins = intlockbins(1:end-1)+0.5*diff(intlockbins([1 2]));
%% Calculate (power normalized) resultant vector
phmag = zeros(numlockbins,numcells);
phangle = zeros(numlockbins,numcells);
for bb = 1:numlockbins
    rvect = cellfun(@(X,Y) nanmean(exp(1i.*X(Y==bb)),1),phase4spike,binintt);
    intlockmag(bb,:) = abs(rvect);
    intlockangle(bb,:) = angle(rvect);
end

%%
% for cc = 1:numcells
% figure
%     plot(inttime4spike{cc},phase4spike{cc},'k.')
%     hold on
%     plot(inttime4spike{cc},phase4spike{cc}+2*pi,'k.')
%     pause
% end
%%
% for cc = 1:numcells
% figure
%     imagesc(inthistbins,phasehistbins,intphasehist(:,:,cc)')
%     hold on
%     imagesc(inthistbins,phasehistbins+2*pi,intphasehist(:,:,cc)')
%     %plot(inthistbincenters,phangle(:,cc),'w')
%     plot(intlockbins,phangle(:,cc)+2*pi,'w')
%     plot(intlockbins,phmag(:,cc).*2*pi-pi,'r')
%     plot(inttime4spike{cc},phase4spike{cc},'w.')
%     plot(inttime4spike{cc},phase4spike{cc}+2*pi,'w.')
%         ylim([-pi 3*pi])
%         axis xy
%     pause
% 
% end

end