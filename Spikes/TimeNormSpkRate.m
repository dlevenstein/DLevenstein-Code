function [meaneventspkrate,t_scaled,alleventspkrate] = TimeNormSpkRate(spiketimes,ints,tbins,restrictint)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% load('Database/BWData/BWRat19_032513/BWRat19_032513_GoodSleepInterval.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_UPSpikeStatsE.mat')
% load('Database/BWData/BWRat19_032513/BWRat19_032513_SSubtypes.mat')
% 
% restrictint = GoodSleepInterval;
% 
% numcells = length(Se);
% numIcells = length(Si);
% for c = 1:numcells
%     CellSpikes{c} = Range(Se{c},'s');
% end
% for c = 1:numIcells
%     ICellSpikes{c} = Range(Si{c},'s');
% end
% 
% ints = [isse.intstarts isse.intends];
% 
% tbins = 100;

%%
numcells = length(spiketimes);
goodtime = [Start(restrictint,'s') End(restrictint,'s')];
intstarts = ints(:,1);
intends = ints(:,2);


numints = length(intstarts);

intstarts = intstarts(intstarts>goodtime(1) & intstarts<=goodtime(2));
intends = intends(intends>goodtime(1) & intends<=goodtime(2));

if intends(1)<intstarts(1)
    intends(1) = [];
end
if intends(end)<intstarts(end)
    intstarts(end) = [];
end


numints = length(intstarts);

numtbins = 3*tbins; %more bins on each side!
t_scaled = linspace(-1,2,numtbins);
%%
display('Scaling Spike Rates');
alleventspkrate = zeros(numtbins,numcells,numints);
for u = 1:numints
    %This is kind of silly/roundabout... interval window with equal window
    %on each size
    ustart = intstarts(u);
    uend = intends(u);
    ulen = uend-ustart;
    ustart = ustart-ulen;
    uend = uend+ulen;
    ulen = uend-ustart;

    dt_u = ulen/numtbins;
    T_u = [0 ustart ulen];
    [spikemat] = SpktToSpkmat(spiketimes,T_u,dt_u);
    while (length(spikemat(:,1))-numtbins)>0
        spikemat(end,:) = [];
    end
    alleventspkrate(:,:,u) = spikemat./dt_u;
end

meaneventspkrate = mean(alleventspkrate,3);


%% Figure
% figure
%     imagesc(t_scaled,1:numcells,log10(meaneventspkrate)')
%     xlim([-0.5 1.5]);

end

