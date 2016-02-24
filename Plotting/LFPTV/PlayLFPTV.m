function [  ] = PlayLFPTV(globalmap,lfpmap,...
    spkidx,stateints,eventints,plotparms,tparms)
%PlayLFPTV (inputs) 
%% DEV


%%
[~,mu,sig] = zscore(log10(globalmap.spec)');

statenames = fieldnames(stateints);
numstates = length(statenames);

eventnames = fieldnames(eventints);
numevents = length(eventnames);

cellgroupnames = fieldnames(spkidx); 
numcellgroups = length(cellgroupnames);

%%
%idea!  loop only over xlim....
%need to downsample LFP 
%% Figure
figure
for t1 = tparms.vidtime(1):tparms.framedt:tparms.vidtime(2)
%t1 = tparms.vidtime(1);

viewwin = [t1 t1+tparms.framewidth];
viewwin_LFP = [round(viewwin(1)*lfpmap.sf):round(viewwin(2)*lfpmap.sf)]+1;
viewwin_LFP = lfpmap.t_LFP>=viewwin(1) & lfpmap.t_LFP<=viewwin(2);

%Global Recording Ref Map
    subplot(12,2,[1,3])
        imagesc(globalmap.t,log2(globalmap.freqs),log10(globalmap.spec))
        axis xy
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        caxis([3.5 6.5])
        caxis([min(mu)-3.5*max(sig) max(mu)+3.5*max(sig)])
        %xlim(viewwin)
        %colorbar('east')
        ylim([log2(globalmap.freqs(1)) log2(globalmap.freqs(end))+0.2])
        %ylabel({'LFP - FFT','f (Hz)'})
    
    hold on
    for ss = 1:numstates
        stateint = stateints.(statenames{ss});
        plot(stateint',log2(globalmap.freqs(end))*ones(size(stateint')),plotparms.statecolors{ss},'LineWidth',8)
    end
        plot([t1 t1],get(gca,'ylim'),'w')
            
        %ylim([0.99 1.01])
        %xlim(vidtime)
        set(gca,'XTickLabel',[]);set(gca,'XTick',[]);
        set(gca,'YTickLabel',[])

%Wavelet Spectrogram    
    subplot(6,1,2:4)
        imagesc(lfpmap.t_LFP(viewwin_LFP),log2(lfpmap.freqs),lfpmap.spec(:,viewwin_LFP))
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        set(gca,'XTickLabel',[])
        xlim(viewwin)
        caxis([-1 4])
        axis xy
        ylabel({'LFP - FFT','f (Hz)'})

%LFP and Spikes     
    subplot(6,1,5:6)    
        hold on
    for ss = 1:numstates
        stateint = stateints.(statenames{ss});
        lineheight = 4;
        plot(stateint',lineheight*ones(size(stateint')),...
            plotparms.statecolors{ss},'LineWidth',3)
    end
    for ee = 1:numevents
        lineheight = 3.5-0.1*ee;
        eventint = eventints.(eventnames{ee});
        plot(eventint',lineheight*ones(size(eventint')),...
            plotparms.eventcolors{ee},'LineWidth',1)
    end
        plot(lfpmap.t_LFP(viewwin_LFP),lfpmap.LFP(viewwin_LFP),'k')

        xlim(viewwin)
    
    hold on
    groupoffset = -35;
    cellnumscale = 10;
    for cc = 1:numcellgroups
        groupspks = spkidx.(cellgroupnames{cc});
        plot(groupspks(:,1),(groupoffset-groupspks(:,2))./cellnumscale,...
            '.','Color',plotparms.spkgroupcolors{cc},'MarkerSize',8)
        groupoffset = groupoffset-max(groupspks(:,2));
    end
        
        xlim(viewwin);ylim([(groupoffset-1)/cellnumscale 4]);
        xlabel('t (s)');
        ylabel('Neuron, Sorted by Rate')
            
        pause(tparms.pausetime)
end

end

