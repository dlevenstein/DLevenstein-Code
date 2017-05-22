function [ output_args ] = bz_LFPSpecToExternalVar( LFP,extvar,specparms,figparms )
%[ output_args ] = bz_LFPSpecToExternalVar( LFP,extvar )
%
%INPUT
%   -LFP        LFP structure in buzcode format:
%               a structure with fields lfp.data, lfp.timestamps
%   -extvar     behavor structure in buzcode format:
%               a structure with fields extvar.data, extvar.timestamps
%               timestamps are in seconds, and aligned
%   -specparms  structure of parameters for the spectrogram
%       .freqs
%
%
%DLevenstein 2017
%%
% LFP = lfp;
% extvar = pupildilation;
% nfreqs = 100;
% specparms.frange = [1 128];
% specparms.nfreqs = 100;
% specparms.ncyc = 5;
% specparms.space = 'log';
% specparms.numvarbins = 10;
% 
% extvar.timestamps = extvar.t_interp;
% extvar.data = extvar.puparea_pxl./max(extvar.puparea_pxl);

%% Calcualte the Wavelet Transform
[freqs,t,spec] = WaveSpec(single(LFP.data),...
    specparms.frange,specparms.nfreqs,specparms.ncyc,...
    1/LFP.samplingRate,specparms.space);
spec = log10(abs(spec));

%% Interpolate time
%Get nearest value of extvar for each point in the LFP
interpvar_LFP = interp1(extvar.timestamps,extvar.data,LFP.timestamps);


%% Binned Spectrogram

binedges = linspace(0,1,specparms.numvarbins+1);
varbins = binedges(1:specparms.numvarbins)+diff(binedges([1 2]));
specbyvar_mean = zeros(specparms.nfreqs,specparms.numvarbins);
specbyvar_std = zeros(specparms.nfreqs,specparms.numvarbins);
for bb = 1:specparms.numvarbins
    inbinidx = interpvar_LFP>=binedges(bb) & interpvar_LFP<=binedges(bb+1);
    inbinspec = spec(:,inbinidx);
    specbyvar_mean(:,bb) = mean(inbinspec,2);
    specbyvar_std(:,bb) = std(inbinspec,[],2);
end

%% correlation for each frequency
sigthresh = 0.05;
[dataLFPcorr,corrp] = corr(interpvar_LFP(~isnan(interpvar_LFP)),spec(:,~isnan(interpvar_LFP))','type','spearman');
corrsig = corrp<=sigthresh;
%%
varcolormap=makeColorMap([0 0 0],[0.8 0 0],specparms.numvarbins);

figure
subplot(2,1,1)
    plot(extvar.timestamps,extvar.data-1,'k')
    hold on
    imagesc(t,log2(freqs),spec)
    axis tight
    %hold on
    axis xy
    LogScale('y',2)
subplot(2,3,4)
    plot(log2(freqs),dataLFPcorr)
    hold on
    plot(log2(freqs(corrsig)),dataLFPcorr(corrsig),'.r')
    plot(log2(freqs([1 end])),[0 0],'r--')
    LogScale('x',2)
    xlabel('Frequency (f)');ylabel('Pupil-Power Correlation (rho)')
    axis tight
subplot(2,3,6)
    imagesc(varbins,log2(freqs),zscore(specbyvar_mean))
    axis xy
    LogScale('y',2)
subplot(2,3,5)
    set(gca,'colororder',varcolormap)
    hold all
    plot(log2(freqs),specbyvar_mean)
    axis tight
NiceSave(['SpecByVar',figparms.plotname],figparms.figfolder,figparms.recname)


end

