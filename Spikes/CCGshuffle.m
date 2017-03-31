function [ ccg,ccgt,ccgnorm,shufflemean,shufflestd ] = CCGshuffle( spiketimes,shufftype,numshuff,CCGparms )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%CCGparms
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%    =========================================================================
%
%DLevenstein Feb2017 - wrapper on FMA toolbox 
%%
%binsize = 0.01;
[ccg,ccgt] = CCG(spiketimes,[],CCGparms{:});

%%
ccg_shuff = zeros([size(ccg) numshuff]);
for ss = 1:numshuff
    %ss
    switch shufftype
        case 'CellID'
            [ spiketimes_shuffle ] = ShuffleSpikelabels(spiketimes);
        case 'jitter'
            jitterwin = 0.15; %...need an input for this... see also
            %jonathan's work showing that spike times should be jittered
            %with in a time bin, not around the center of a spike
             [ spiketimes_shuffle ] = JitterSpiketimes(spiketimes,jitterwin);
    end
    ccg_shuff(:,:,:,ss) = CCG(spiketimes_shuffle,[],CCGparms{:});
end

shufflemean = mean(ccg_shuff,4);
shufflestd = std(ccg_shuff,[],4);

ccgnorm = (ccg-shufflemean)./shufflestd;
ccgsig = abs(ccgnorm)>3;
%% Figure for example CCG
% numcells = length(spiketimes);
% figure
%     excell = randi(numcells);
% excell2 = randi(numcells);
% subplot(6,2,1)
%     bar(ccgt,squeeze(ccg(:,excell,excell2)),'facecolor','k')
%     hold on
%     plot(ccgt,squeeze(shufflemean(:,excell,excell2)),'r-','linewidth',2)
%     plot(ccgt,squeeze(shufflemean(:,excell,excell2)+shufflestd(:,excell,excell2)),'r--','linewidth',2)
%     plot(ccgt,squeeze(shufflemean(:,excell,excell2)-shufflestd(:,excell,excell2)),'r--','linewidth',2)
%     axis tight
%     box off
%     plot([0 0],get(gca,'ylim'),'k--')
%    %  xlim([-0.2 0.2])
% subplot(6,2,[3 5])
%     %colormap(gca,'jet')
%     %[~,refsort] = sort(CCGpeakdelay(excell,:));
%     imagesc(ccgt,[1 numcells],(squeeze(ccgnorm(:,excell,:)))')
%     hold on
%     plot([0 0],[0 numcells+1],'w--')
%    % plot(get(gca,'xlim'),find(refsort==excell).*[1 1],'k')
%    % ColorbarWithAxis([-5 5],'Z Score')
%     
