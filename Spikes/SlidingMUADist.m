function [MUAhist,t_hist,histbins,MUAhist_sm,sm_histbins,MUA,t_spkmat] = SlidingMUADist(cellspikes,histparms,T)
%[MUAhist,t_hist,histbins,MUAhist_sm,sm_histbins] = SlidingMUADist(cellspikes,histparms,T)
%
%INPUTS
%   cellspikes:     cell array of spike times for all cells
%   histparms:      struct of time/histogram parms
%       .dt         time binning for histogram output (s)
%       .winlength  time window to count over for histogram (s)
%       .spkbinsize time bin in which to count spikes (s)
%       .smoothwin  time to smooth MUA over
%       .histsmbins factor to increase smooth histogram by
%       .numbins
%   T:              time window (passed to SpktToSpkmat)
%                   T = [t_start t_offset t_end]
%                   t_offset: time t = 0
%                   t_start: start time relative to t_offset
%                   t_end: end time relative to t_offset
%
%
%
%TO DO
%   -Add option for MUA (if one "cell" of spike times, normalize FR between
%   min and max instead of number of cells firing)
%   -might want to make it so raw hist is not smoothed in time as well...
%   but MUAhist_sm is?
%   -Boxcar for MUA...
%
%Updated 9/27/15.... did a lot of sloppy stuff.  Sorry future me :(
%DLevenstein

%% Parms

% histparms.dt = 0.5;
% histparms.winlength = 15; %s
% histparms.spkbinsize = 0.01; %
% histparms.smoothwin = 0.1;
% histparms.histsmbins = 4;


dt_hist = histparms.dt_hist;
dt_MUA = histparms.dt_MUA;
winlength = histparms.winlength;
spkbinsize = histparms.spkbinsize;
smoothwin = histparms.smoothwin;
histsmbins = histparms.histsmbins;

tfact= round(dt_hist/dt_MUA);

numcells = length(cellspikes);

numwinbins = ceil(winlength/dt_MUA);
span_smoothMUA = ceil(smoothwin/dt_MUA);
%%
% load('Database/BWData/BWRat19_032513/BWRat19_032513_SSubtypes.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_GoodSleepInterval.mat')
% clear cellspikes 
% clear ICellSpikes
% 
% numcells = length(Se);
% numIcells = length(Si);
% for c = 1:numcells
%     cellspikes{c} = Range(Se{c},'s');
% end
% for c = 1:numIcells
%     ICellSpikes{c} = Range(Si{c},'s');
% end
% cellspikes = [cellspikes,ICellSpikes];
% numcells = numIcells+numcells;
%%

%T = [0 End(GoodSleepInterval,'s')];
overlap = spkbinsize/dt_MUA; 
[spikemat,t_spkmat,spindices] = SpktToSpkmat(cellspikes, T, dt_MUA,overlap);
t_spkmat = t_spkmat-0.5*dt_MUA;


%Don't double the same cell spiking in a window
%FOR PERC CELLS
%spikemat = spikemat>0;
MUA = sum(spikemat,2)/numcells;
%FOR FR



MUA = smooth(MUA,span_smoothMUA);

histwins = repmat(MUA,1,numwinbins);
for w = 1:numwinbins
    histwins(:,w) = circshift(histwins(:,w),-(w-1)+0.5*numwinbins);
end


% Ehistbins = linspace(0,max(EMUA),numhistbins);

if isfield(histparms,'numbins')
    rawhistbins = histparms.numbins;
else
    rawhistbins = (numcells+1);
end
histbins = linspace(0,1,rawhistbins);

%histbins = linspace(0,max(MUA),20);

sm_histbins = linspace(0,1,histsmbins*rawhistbins);
% nhistbins = 1000;
% binscale = nhistbins/rawhistbins;

t_hist = t_spkmat(1:tfact:end);
MUAhist = hist(histwins(1:tfact:end,:)',histbins)'./numwinbins;
MUAhist_sm = hist(histwins(1:tfact:end,:)',sm_histbins)'*histsmbins;
for i_t = 1:length(t_hist)
    MUAhist_sm(i_t,:) = smooth(MUAhist_sm(i_t,:),2*histsmbins);
% 
%     tol = 0.0001;
%     newhist = [1 1];
%     while length(newhist) ~= nhistbins
%         [P,Q] = rat(binscale,tol);
%         newhist = resample(EMUAhist(i_t,:),P,Q,'linear');
%         tol = tol/10;
%     end
%     EMUAhist(i_t,:) = newhist;
end


% if isfield(histparms,'numbins')
%     MUAhist = MUAhist*numcells;
% else

 %% Figure
% viewwin = [378 2300];
% %viewwin = [1200 1600];
% viewwin = [1375 1450];
% %viewwin=[3000 3008];
% 
% figure
%     subplot(4,1,1)
%         plot(spindices(:,1),spindices(:,2),'.')
%         xlim(viewwin);ylim([-numIcells numcells+1]);
%         xlabel('t (s)');
%         ylabel('Neuron')
%     subplot(4,1,2)
%         plot(t_spkmat,MUA,'k-')
%         xlim(viewwin)
%     subplot(4,1,3)
%         imagesc(t_hist,histbins,MUAhist')
%         xlim(viewwin)
%         %caxis([0 1])
%         ylim([histbins(1)-0.5*histbins(2) max(MUA)])
%         %colorbar
%         axis xy
%     subplot(4,1,4)
%         imagesc(t_hist,sm_histbins,MUAhist_sm')
%         xlim(viewwin)
%         %caxis([0 1])
%         ylim([sm_histbins(1)-0.5*sm_histbins(2) max(MUA)])
%         %colorbar
%         axis xy
% 
%          colormap('pink')
% 
%          
 %% Figure: Example Histograms
% viewwin = [378 2300];
% ex1win = [700 715];
% ex2win = [1390 1405];
% ex3win = [1185 1200];
% 
% figure
% 	subplot(4,1,1)
%         imagesc(t_hist,sm_histbins,MUAhist_sm')
%         xlim(viewwin)
%         %caxis([0 1])
%         ylim([histbins(1) max(MUA)])
%         %colorbar
%         axis xy
% 
%          colormap('pink')
%           
% 	subplot(4,3,[4:5])
%         plot(t_spkmat,MUA,'k-')
%         xlim(ex1win)
% 	subplot(4,3,6)
%         hold on
%         bar(histbins,MUAhist(find(t_hist>ex1win(1),1),:))
%         plot(sm_histbins,MUAhist_sm(find(t_hist>ex1win(1),1),:))
%         xlim([0 0.1])
%         
% 	subplot(4,3,[7:8])
%         plot(t_spkmat,MUA,'k-')
%         xlim(ex2win)
% 	subplot(4,3,9)
%         hold on
%         bar(histbins,MUAhist(find(t_hist>ex2win(1),1),:))
%         plot(sm_histbins,MUAhist_sm(find(t_hist>ex2win(1),1),:))
%         xlim([0 0.1])
%         
% 	subplot(4,3,[10:11])
%         plot(t_spkmat,MUA,'k-')
%         xlim(ex3win)
% 	subplot(4,3,12)
%         hold on
%         bar(histbins,MUAhist(find(t_hist>ex3win(1),1),:))
%         plot(sm_histbins,MUAhist_sm(find(t_hist>ex3win(1),1),:))
%         xlim([0 0.1])
        
end

