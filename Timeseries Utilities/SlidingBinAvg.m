function [binavg,bin_t] = SlidingBinAvg(magnitudes,times,dt,overlap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
% magnitudes = dnupmag;
% times = dnuptimes;
% dt = 1;
% overlap = 6;

winsize = dt*overlap;

maxt = max(times);
numbins = ceil(maxt/dt);
maxt = numbins*dt;
binedges = linspace(0,maxt,numbins+1);


[bin_t,binsums,~,binnum] = BinDataTimes(magnitudes,times,binedges,'sum');
%[numeventsmat,tbins] = SpktToSpkmat({times},maxt,dt,overlap);

%Should probably check that FConv has the convolution centered correctly,
%might be off by a bin...
binsumswin = FConv(ones(overlap,1),binsums);
binnumwin = FConv(ones(overlap,1),binnum);
binsumswin(binnumwin<0.5)=0;
binnumwin(binnumwin<0.5)=0;
%numeventswin = FConv(ones(overlap,1),numeventsmat);
%binsums = repmat(binsums,1,overlap);

binavg = binsumswin./binnumwin;

%%
% figure
%     hold on
%     plot(times,magnitudes,'.')
%     plot(bincenters,binnumwin)
%     plot(bincenters,binsumswin)
%     plot(bincenters,binavg)

end

