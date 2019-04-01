function [ amp1bins,phasebins,sig2powerskew,sig2prefangle,phaseamphist,...
    binsig,threshsig] = ...
    PhaseAmpCouplingByAmp(sig1phase,sig1amp,sig2amp,varargin)
%PhaseAmpCouplingByAmp(sig1phase,sig1amp,sig2amp,ampbins,numphasebins)
%calculates phase amplitude coupling between the phase of signal 1 and the
%amplitude of signal 2, with respect to the amplitude of signal 1.
%
%INPUTS
%   sig1phase
%   sig1amp
%   sig2amp
%
%   (options)
%       'numAmpbins'    default: 10
%       'numPhbins'     default: 20
%       'AmpBounds'     default: [-2.5 2.5] for Z-normed amp
%       'AmpZNorm'      default: true
%       'showFig'       default: false
%
%OUTPUTS
%   binsig          bin/threshsig are two structures with coupling and 
%   threshsig       significance for each power bin. Threshsig has _above
%                   and _below for many fields.
%       .powerbin       power bins
%       .coupling       coupling strength
%       .perc           coupling value of 99.5% shuffle (above this is "significant")
%       .sigstdbin      number of standard deviations coupling is from mean of shuffles 
%       .sigthresh      power threshold aabove which there is significant coupling
%       .sig            t/f is it significant?
%
%DLevenstein Fall 2016
%%
p = inputParser;

addParameter(p,'numAmpbins',10)
addParameter(p,'numPhbins',20)
addParameter(p,'AmpBounds',[-2.5 2.5])
addParameter(p,'AmpZNorm',true)
addParameter(p,'shufflesig',true)
addParameter(p,'showFig',false)

parse(p,varargin{:})
numAmpbins = p.Results.numAmpbins;
numPhbins = p.Results.numPhbins;
AmpBounds = p.Results.AmpBounds;
AmpZNorm = p.Results.AmpZNorm;
SHOWFIG = p.Results.showFig;
shufflesig = p.Results.shufflesig;
threshsig = p.Results.threshsig;


%%
phasebins = linspace(-pi,pi,numPhbins+1);
phasebins=phasebins(1:end-1)+diff(phasebins(1:2));
amp1bins = linspace(AmpBounds(1),AmpBounds(2),numAmpbins);
amp2bins = linspace(0,5,numAmpbins);

%%
if AmpZNorm
    sig1amp = zscore(sig1amp);
end
sig2amp = sig2amp./mean(sig2amp);
%%
sig1binpower = interp1(amp1bins,amp1bins,sig1amp,'nearest');
sig1binphase = interp1(phasebins,phasebins,sig1phase,'nearest');


powerhist = zeros(numAmpbins,numPhbins);
sig2prefangle = zeros(numAmpbins,1);
sig2powerskew = zeros(numAmpbins,1);
for bb = 1:numAmpbins
    thisampbin = sig1binpower==amp1bins(bb);
    for bbb = 1:numPhbins
        ampwintimes = sig2amp(thisampbin & sig1binphase==phasebins(bbb));
        phaseamphist(bb,bbb) = mean(log10(ampwintimes));
    end
    sig2powerskew(bb) = mean(sig2amp(thisampbin).*exp(1i.*sig1phase(thisampbin)));
    sig2prefangle(bb) = angle(sig2powerskew(bb));
    sig2powerskew(bb) = abs(sig2powerskew(bb));
    
    %Threshold test
	abovethresh = sig1amp>=amp1bins(bb);
	belowthresh = sig1amp<=amp1bins(bb);
	abovethreshcoupling(bb) = abs(mean(sig2amp(abovethresh).*exp(1i.*sig1phase(abovethresh))));
	belowthreshcoupling(bb) = abs(mean(sig2amp(belowthresh).*exp(1i.*sig1phase(belowthresh))));
    
    if shufflesig
        numshuff = 200;
        for ss = 1:numshuff
            shuffamps = randsample(sig2amp(thisampbin),length(sig2amp(thisampbin)),true);
            shuffskew(ss) = mean(shuffamps.*exp(1i.*sig1phase(thisampbin)));
            shuffskew(ss) = abs(shuffskew(ss));
            
            shuffamps = randsample(sig2amp(abovethresh),length(sig2amp(abovethresh)),true);
            shuffabove(ss) = abs(mean(shuffamps.*exp(1i.*sig1phase(abovethresh))));
            
            shuffamps = randsample(sig2amp(belowthresh),length(sig2amp(belowthresh)),true);
            shuffbelow(ss) = abs(mean(shuffamps.*exp(1i.*sig1phase(belowthresh))));

        end
            shuffmean(bb) = mean(shuffskew);
            shuffstd(bb) = std(shuffskew);
            

            perc_bin(bb) = max(shuffskew);            
            perc_above(bb) = max(shuffabove);
            perc_below(bb) = max(shuffbelow);
            
            sigstdbin(bb) = (sig2powerskew(bb)-shuffmean(bb))./shuffstd(bb);
            sigstdabove(bb) = (abovethreshcoupling(bb)-mean(shuffabove))./std(shuffabove);
            sigstdbelow(bb) = (belowthreshcoupling(bb)-mean(shuffbelow))./std(shuffbelow);
    end
    
  
end

sig_above = abovethreshcoupling>perc_above;
sig_below = belowthreshcoupling>perc_below;
sig_bin = sig2powerskew>perc';

sigthresh_type2 = max(amp1bins(sig_above & ~sig_below));
sigthresh_type1 = min(amp1bins(sig_bin));
%%
binsig.powerbin = amp1bins;
binsig.coupling = sig2powerskew;
binsig.perc_bin = perc_bin;
binsig.sigstdbin = sigstdbin;
binsig.sigthresh = sigthresh_type1;
binsig.sig_bin = sig_bin;

threshsig.powerbin = amp1bins;
threshsig.coupling_above = abovethreshcoupling;
threshsig.coupling_below = belowthreshcoupling;
threshsig.perc_above = perc_above;
threshsig.perc_below = perc_below;
threshsig.sigstd_above = sigstdabove;
threshsig.sigstd_below = sigstdbelow;
threshsig.sig_above = sig_above;
threshsig.sig_below = sig_below;
threshsig.sigthresh = sigthresh_type2;





%%
% figure
% hist(shuffskew)
%% Figure
if SHOWFIG
    rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
    plotx = linspace(-pi,3*pi,100);
    figure
    subplot(2,2,1)
    hold on
        imagesc(phasebins,amp1bins,phaseamphist)
        imagesc(phasebins+2*pi,amp1bins,phaseamphist)
        plot(sig2prefangle,amp1bins,'.k')
        plot(sig2prefangle+2*pi,amp1bins,'.k')
        plot(plotx,bz_NormToRange(cos(plotx),AmpBounds),'k')
        colormap(gca,rwbcolormap)
        axis xy
        axis tight
        ColorbarWithAxis([-1 1],['Mean Amp.'])
        %caxis([0 2])
      %  xlim([-pi 3*pi]);ylim(amp1bins([1 end]))
        xlabel('Signal 1 Phase');ylabel('Signal 1 Amp (Z)')
    subplot(2,4,3)
        plot(sig2powerskew,amp1bins,'k','LineWidth',1)
        hold on;box off
        plot(perc,amp1bins,'r--','LineWidth',1)
        ylabel('Signal 1 Amp. (Z)');
        xlabel('Phase-Amp. Modulation (mrl)')
        axis tight
    subplot(2,4,4)
        plot(abovethreshcoupling,amp1bins,'k','LineWidth',1)
        %plot(sigabove,amp1bins,'k','LineWidth',1)
        hold on;box off
        %plot(sigbelow,amp1bins,'r','LineWidth',1)
        plot(belowthreshcoupling,amp1bins,'r','LineWidth',1)
        plot(perc_above,amp1bins,'k--','LineWidth',1)
        plot(perc_below,amp1bins,'r--','LineWidth',1)
        ylabel('Signal 1 Amp. (Z)');
        xlabel('Phase-Amp. Modulation (mrl)')
        axis tight
    subplot(4,2,6)
        histogram(sig1amp,amp1bins)
        xlabel('Signal 1 Amp. (Z)');
        ylabel('Occupancy')
        axis tight
        title('Signal1/2 Amp. Distributions')

    subplot(4,2,8)
        histogram(sig2amp,amp2bins)
        xlabel('Signal 2 Amp. (Z)');
        ylabel('Occupancy')
        axis tight 

    subplot(2,2,3)
        gasphist = hist3([sig1amp,sig2amp],{amp1bins,amp2bins});
        imagesc(amp1bins,amp2bins,gasphist')
        axis xy
        xlabel('Signal 1 Power');ylabel('Signal 2 Power')
end



end
