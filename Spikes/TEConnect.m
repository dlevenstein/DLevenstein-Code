function [peakTE,CI,sigTE,allTE,maxdelay] = TEConnect(spikes,int,timerange,SIGTEST)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%See: Shimono and Beggs 2014, Ito et al 2011
%https://code.google.com/p/transfer-entropy-toolbox/
%
%%
SELFCONNECT = false;
%SIGTEST = true;
numshuff = 100;

if ~exist('timerange','var')
    timerange = [1 100];
end 
delays = [timerange(1):timerange(2)];


    
%Convert Spiketimes to a spike matrix with 1ms bins
dt = 0.001;
[spikemat,t] = SpktToSpkmat(spikes, [], dt);

%Restrict to spikes in int (if applicable)
if exist('int','var')
    intidx = INTtoIDX(int,length(t),1/dt);
    intidx = logical(intidx);
    spikemat = spikemat(intidx,:);
end

%Convert Within-State Spikemat to asdf for TE function
asdf = SparseToASDF(spikemat',dt*1000);

%Calculate Transfer Entropy
i_order =3; %postsynaptic
j_order =3; %presynaptic
i_order =1; %postsynaptic
j_order =1; %presynaptic

[peakTE,CI,allTE,maxdelay] = ASDFTE_parallel(asdf, delays,i_order, j_order);

%%
if ~SELFCONNECT
    peakTE(1:size(peakTE,1)+1:end) = nan;
    CI(1:size(CI,1)+1:end) = nan;
end

shuffTE = zeros([size(peakTE),numshuff]);
%shuffCI = zeros([size(CI),numshuff]);
if SIGTEST
    parfor tt = 1:numshuff
        if mod(tt,5)==1
            %display(['Shuffle ',num2str(tt),' of ',num2str(numshuff)])
        end
        [shuffTE(:,:,tt)] = TEConnect(JitterSpiketimes(spikes,0.02),int,timerange,false);
    end
end

%%
TEz = zscore(cat(3,peakTE,shuffTE),[],3);
%TEz = TEz(:,:,1);
sigTE = TEz(:,:,1)>3;

%%
if SIGTEST
figure
hold on
    plot(log10(peakTE(:)),CI(:),'k.','MarkerSize',15)
    plot(log10(peakTE(sigTE)),CI(sigTE),'r.','MarkerSize',15)
    %plot(log10(shuffTE(:)),shuffCI(:),'r.','MarkerSize',10)
    xlabel('Transfer Entropy (bits)')
    LogScale('x',10)
    ylabel('Confidence Index')
    legend('All Pairs','Pair with TE > 3SD shuffled')
end
%%
% figure
% hold on
%     plot(log10(peakTE(28,2)),CI(28,2),'k.','MarkerSize',15)
%     plot(log10(squeeze(shuffTE(28,2,:))),squeeze(shuffCI(28,2,:)),'r.','MarkerSize',15)
%     %plot(log10(shuffTE(:)),shuffCI(:),'r.','MarkerSize',10)
%     xlabel('Transfer Entropy (bits)')
%     LogScale('x',10)
%     ylabel('Confidence Index')
%     %legend('All Pairs','Pair with TE > 3SD shuffled')
end

