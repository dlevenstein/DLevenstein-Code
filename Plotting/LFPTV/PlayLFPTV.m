function [ output_args ] = PlayLFPTV( input_args )
%PlayLFPTV (inputs) 
%% DEV


%%
%% Figure
figure
for t1 = vidtime(1):framesize:vidtime(2)


%t1 = 0;
viewwin = [t1 t1+windowsize];
viewwin_LFP = [round(viewwin(1)*sf_LFP):round(viewwin(2)*sf_LFP)]+1;
viewwin_LFP = t_LFP>=viewwin(1) & t_LFP<=viewwin(2);

    subplot(18,1,1)

        plot(REM',ones(size(REM')),'g',...
            SWS',ones(size(SWS')),'b',...
            Wake',ones(size(Wake')),'r',...
            t1,1,'w.',...
            'LineWidth',10)
        ylim([0.99 1.01])
        xlim(vidtime)
        set(gca,'XTickLabel',[]);set(gca,'XTick',[]);
        set(gca,'YTickLabel',[])

    subplot(16,1,10:12)
        plot(REM', 2.4*ones(size(REM')),'g',...
            SWS', 2.4*ones(size(SWS')),'b',...
            Wake', 2.4*ones(size(Wake')),'r',...
            UPs(:,1),2.3*ones(size(UPs(:,1))),'+g',...
            UPs(:,2),2.3*ones(size(UPs(:,2))),'+r',...
            SPs(:,1),2*ones(size(SPs(:,1))),'+c',...
            SPs(:,2),2*ones(size(SPs(:,2))),'+m',...
            t_LFP(viewwin_LFP),LFP(viewwin_LFP),'k')

        xlim(viewwin)
        ylim([-3 3])

    
    subplot(16,1,2:9)
        imagesc(t_LFP(viewwin_LFP),log2(freqs),zspec(:,viewwin_LFP))
        set(gca,'YTick',(log2([1 2 4 8 16 32 64 128])))
        set(gca,'YTickLabel',{'1','2','4','8','16','32','64','128'})
        set(gca,'XTickLabel',[])
        xlim(viewwin)
        caxis([-1 4])
        axis xy

        
    subplot(4,1,4)
        plot(spindices(:,1),spindices(:,2),'.',...
            Ispindices(:,1),-Ispindices(:,2),'r.','MarkerSize',6)
        xlim(viewwin);ylim([-numIcells numcells+1]);
        xlabel('t (s)');
        ylabel('Neuron, Sorted by Rate')
        
        
        pause(pausetime)
end

end

