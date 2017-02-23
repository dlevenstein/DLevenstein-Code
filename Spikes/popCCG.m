function [popccg,ccgt] = popCCG( spiketimes,CCGparms )
%[popccg,ccgt] = popCCG( spiketimes,CCGparms )
%
%CCGparms
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%    =========================================================================
%
%OUTPUT
%   popccg  population CCG in units of spks/spk/cell/s
%
%DLevenstein Feb2017 - wrapper on FMA toolbox 
%%

%Calculate all ccgs
[allccg,ccgt] = CCG(spiketimes,[],CCGparms{:});




%%
numcells = length(spiketimes);
binsize = diff(ccgt(1:2));
popccg = zeros(length(ccgt),numcells);

for cc = 1:numcells
    %Remove ACG
    allccg(:,cc,cc) = zeros(size(allccg(:,cc,cc)));
    %Add the CCGs of all other cells around reference cell
    popccg(:,cc) = sum(allccg(:,cc,:),3);
    
    %Normalize to number of spikes,cells,timebin
    numspks = length(spiketimes{cc});
    popccg(:,cc) = popccg(:,cc)./numspks./(numcells-1)./binsize;
end

%%
figure
imagesc(ccgt,[0 1],(popccg(:,sortccgCM))')
hold on
plot([0 0],[0 1],'w--')
plot(sort(ccgCM),(1:numcells)/numcells,'r')
axis xy
colorbar
xlim([-0.2 0.2])

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
