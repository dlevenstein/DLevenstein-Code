function [Spindles,ChannelProperties] = DetectSPINDLES(lfpfilename,channels2use,NREMints,detectionchannel,figfolder,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%varargin: 'detectionchannel','phasereferencechannel'. default is to select
%this...
%otherparms
%sf  (default = 1250Hz)
%channels2use (0-indexed as in neuroscope)
%
%outputs to add: phasereferencechannel, spindlephase
%
%
%INPUT
%   lfpfilename     filepath of the .lfp/.eeg file
%   channels2use    list of (good) cortical channels on the same shank,
%                   same 0-indexing as in neuroscope
%   NREMints        [Nints x 2] start and end times of NREM time periods
%   detectchannel   'best' - selects best channel for detection
%                   'all' detects spindles on all channels
%                   -or can be vector of channel numbers, subset of the
%                   channels2use
%   figfolder       folder for output figures
%                   
%
%
%OUTPUT
%   Spindles
%       .ints       [Nints x 2] start and end times of ID'd spindle
%                   activity
%       .Detectedon {Nints} cell array listing the channels that each
%                   spindle was identified on
%       .detectchannelmat   [Nints x Nchannels] matrix with the same info
%   ChannelProperties
%       .powerskew
%       .powerfreq
%       .channum
%       .gammapower
%       .spindlephase
%% DEV

% lfpfilename = '/Users/dlevenstein/Dropbox/Research/Datasets/AYAData/AYA2_150216/AYA2_150216.eeg';
% channels2use = [225,227,228,229,230,231,232,233,234,235,236,237,238,239,...
%      240,241,242,243,244,245,246,247,248,249,250,251,252,254,255]';
% detectionchannel = 'best';
% phasereferencechannel='detectionchannel';
% figfolder = '/Users/dlevenstein/Dropbox/Research/Current Projects/TransitionsAndSpindles/AnalysisScripts/AnalysisFigs/SpindleDetectAnalysis/';

%% Load the xml for metadata
[filepath,filename] = fileparts(lfpfilename);
xmlfilename = fullfile(filepath,[filename,'.xml']);

%Load Metadata from xml
Par = LoadPar_SleepScore(xmlfilename);
Fs = Par.lfpSampleRate; % Hz, LFP sampling rate


numchannels = length(channels2use);

%% DEV
%Load State Scoring
%load(fullfile(filepath,[filename,'_SleepScore.mat']))
%NREMints = StateIntervals.NREMpacket;

%% Load the LFP for Power Skew Analysis
downsamplefactor = 10;
sf_down = Fs/downsamplefactor;
LFP = LoadBinary_DL(lfpfilename,'frequency',Fs,'nchannels',Par.nChannels,...
    'channels',(channels2use)+1,'duration',max(NREMints(:)),'downsample',downsamplefactor);
t_LFP = [1:length(LFP(:,1))]'./sf_down;

NREMLFP = IsolateEpochs2(LFP,NREMints,0,sf_down);
t_NREMLFP = IsolateEpochs2(t_LFP,NREMints,0,sf_down);
NREMLFP = cat(1,NREMLFP{:});
t_NREMLFP = cat(1,t_NREMLFP{:});
[NREMLFP,NREM_meanLFP,NREM_sigLFP] = zscore(NREMLFP);

%% Correlation between LFP channels
LFPcorr = corr(NREMLFP,'type','spearman');

%% Range for Power Distribution
frange = [6 32];
nfreqs = 50;

meanspec = zeros(nfreqs,numchannels);
medianspec = zeros(nfreqs,numchannels);

%% Calculate Spectrogram for all channels
for cc = 1:numchannels
    
    freqlist = logspace(log10(frange(1)),log10(frange(2)),nfreqs);
    window = 1;
    noverlap = 0.8;
    window = window*sf_down;
    noverlap = noverlap*sf_down;
    [spec,freqs,t_FFT] = spectrogram(NREMLFP(:,cc),window,noverlap,freqlist,sf_down);
    spec = abs(spec);

    meanspec(:,cc) = mean(abs(spec),2);
    medianspec(:,cc) = median(abs(spec),2);
end

%% Spindle Frequency Skew
sprange = [10 17];
meanskew = meanspec./medianspec;
[maxskew,skewfreq] = max(meanskew(freqs>=sprange(1) & freqs<=sprange(2),:));
spfreqs = freqs((freqs>=sprange(1) & freqs<=sprange(2)));
skewfreq = spfreqs(skewfreq);

    ChannelProperties.powerskew = meanskew;
    ChannelProperties.powerfreq = freqs;
    ChannelProperties.channum = channels2use;
    ChannelProperties.maxskew = maxskew;
    ChannelProperties.skewfreq = skewfreq;

% Best channel has most skew in power of filtered signal.
[spLFP,spPOWER,spPHASE] = FiltNPhase(NREMLFP,sprange,sf_down);
spindlebandmean = mean(spPOWER,1);
spindlebandmedian = median(spPOWER,1);
spindlebandskew = spindlebandmean./spindlebandmedian;

[~,bestskewchannelidx] = max(spindlebandskew(1:end-1));
bestskewchannelnumber = channels2use(bestskewchannelidx);

     ChannelProperties.bestskewchannelnumber = bestskewchannelnumber;


%Spindle Power Histogram All Channels
[spowerstats.allhist,spowerstats.powerbins] = hist(spPOWER,100);
spowerstats.allhist = bsxfun(@(x,y) x./y,spowerstats.allhist,max(spowerstats.allhist));

%Spindle Power Correlation All Channels, Phase Coherence
spPowerCorr = corr(spPOWER,'type','spearman');
%spPhaseCoh = 

switch detectionchannel
    case 'best'
        detectchannelidx = bestskewchannelidx;
        detectchannelnumber = bestskewchannelnumber;
    case 'all'
        detectchannelidx = 1:length(channels2use);
        detectchannelnumber = channels2use;
    otherwise
        detectchannelidx = find(channels2use==detectionchannel);
        detectchannelnumber = detectionchannel;
end

numdetectchannels = length(detectchannelnumber);

for cc = 1:numdetectchannels
    currentchannelidx = detectchannelidx(cc);
    currentchannelnumber = detectchannelnumber(cc);
    display(['Detecting Spindles on Channel ',num2str(currentchannelnumber),' (',num2str(cc),' of ',num2str(numdetectchannels),')'])


    %Power Density  and Power-Power Correlation for Best Channel
    [spec,freqs,t_FFT] = spectrogram(NREMLFP(:,currentchannelidx),window,noverlap,freqlist,sf_down);
    spec = abs(spec);

    numpowerbins = 200;
    specpowerbins = linspace(-0.5,1.75,numpowerbins);
    [powerdist_best] = hist(log10(abs(spec))',specpowerbins);
    [colmax,peakf] = max(powerdist_best,[],1);
    powerdist_best = bsxfun(@(X,Y) X./Y,powerdist_best,colmax);

    meanskew_best = meanskew(:,currentchannelidx);
    skewfreq_best = skewfreq(currentchannelidx);
    maxskew_best = maxskew(currentchannelidx);

    %% Spindle Phase-Gamma Power Coupling on best channel
    downsample_best = 2;
    sf_best = Fs./downsample_best;
    bestLFP = LoadBinary_DL(lfpfilename,'frequency',Fs,'nchannels',Par.nChannels,...
        'channels',currentchannelnumber+1,'duration',max(NREMints(:)),'downsample',downsample_best);

    t_LFP_best = [1:length(bestLFP(:,1))]'./sf_best;

    NREMLFP_detectchan = IsolateEpochs2(bestLFP,NREMints,0,sf_best);
    t_NREMLFP_best = IsolateEpochs2(t_LFP_best,NREMints,0,sf_best);
    NREMLFP_detectchan = cat(1,NREMLFP_detectchan{:});
    t_NREMLFP_best = cat(1,t_NREMLFP_best{:});
    NREMLFP_detectchan = zscore(NREMLFP_detectchan);

    %% Filter: Gamma and Spindle

    garange = [80 160];
    [spLFP,spPOWER,spPHASE] = FiltNPhase(NREMLFP_detectchan,sprange,sf_best);
    [gaLFP,gaPOWER,gaPHASE] = FiltNPhase(NREMLFP_detectchan,garange,sf_best);
    
    %Calcluate local spindle-gamma PAC
    spindlegammaPAC = mean(gaPOWER.*exp(1i.*spPHASE)./mean(gaPOWER));
    
    %Z-log transform for normality
    spPOWER = zscore(log10(spPOWER));
    gaPOWER = zscore(log10(gaPOWER));
    %%
    numbins = 21;
    phasebins = linspace(-pi,pi,numbins+1);phasebins=phasebins(1:end-1)+diff(phasebins(1:2));
    powerbins = linspace(-2,2,numbins);

    spbinpower = interp1(powerbins,powerbins,spPOWER,'nearest');
    spbinphase = interp1(phasebins,phasebins,spPHASE,'nearest');

    gammahist = zeros(numbins,numbins);
    gammaangle = zeros(numbins,1);
    gammaskew = zeros(numbins,1);
    for bb = 1:numbins
        for bbb = 1:numbins
            gamatimes = gaPOWER(spbinpower==powerbins(bb) & spbinphase==phasebins(bbb));
            gammahist(bb,bbb) = mean(gamatimes);
        end
        gammaskew(bb) = mean(gaPOWER(spbinpower==powerbins(bb)).*exp(1i.*spPHASE(spbinpower==powerbins(bb))));
        gammaangle(bb) = angle(gammaskew(bb));
        gammaskew(bb) = abs(gammaskew(bb));
    end



    %% Figure: LFP on all channels: Spindle Power 
    winsize = 3; %s
    randwin = t_NREMLFP(randi(length(NREMLFP)))+[0 winsize];
    rwbcolormap = [makeColorMap([0 0 1],[0.8 0.8 1],[1 1 1]);makeColorMap([1 1 1],[1 0.8 0.8],[1 0 0])];
    rwbcolormap2 = makeColorMap([0 0 1],[1 1 1],[1 0 0]);
    rainbowmap = RainbowColors(numchannels);
    distmap = makeColorMap([1 1 1],[0 0.6 0],[0.8 0.7 0]);


    numrepchannels = 15;
    repchannels = round(linspace(1,numchannels,numrepchannels));
    figure
        subplot(3,2,[1 3])
            plot(log2(freqs),bsxfun(@(X,Y) X+Y,(meanskew(:,repchannels)-1).*15,-repchannels),'k','LineWidth',1)
            hold on
            plot(log2(freqs),bsxfun(@(X,Y) X+Y,zeros(size((meanskew(:,repchannels)))),-repchannels+(meanskew(1,repchannels)-1).*15),'--k')
            plot(log2(freqs),bsxfun(@(X,Y) X+Y,(meanskew(:,currentchannelidx)-1).*15,-currentchannelidx),'r','LineWidth',2)
            xlabel('f (Hz)');ylabel('Channel')
            LogScale('x',2)
            axis tight
            title('Power Distirbution Skew by Depth')
            plot(spindlebandskew(1:end)-1+log2(freqs(end)),-[0:numchannels-1],'ro-')
            plot(spindlebandskew(currentchannelidx)-1+log2(freqs(end)),-[currentchannelidx-1],'r*')
            plot(log2(freqs(end)).*[1 1],get(gca,'ylim'),'k')
            plot(log2(sprange(1)).*[1 1],get(gca,'ylim'),'k--')
            plot(log2(sprange(2)).*[1 1],get(gca,'ylim'),'k--')

        subplot(3,3,7)
        colormap(gca,rwbcolormap)
            imagesc(spPowerCorr)
            hold on
            plot(currentchannelidx,currentchannelidx,'k*')
            xlabel('Channel');ylabel('Channel')
            set(gca,'Ytick',[]);set(gca,'XTick',[]);
            caxis([-1 1])
            colorbar
            title('Correlation: Spindle Power')

        subplot(3,3,8)
        colormap(gca,rwbcolormap)
            imagesc(LFPcorr)
            hold on
            plot(currentchannelidx,currentchannelidx,'k*')
            xlabel('Channel');ylabel('Channel')
            set(gca,'Ytick',[]);set(gca,'XTick',[]);
            caxis([-1 1])
            colorbar
            title('Correlation: LFP')

          subplot(3,2,2)
         colormap(gca,distmap)
            plot(log2(freqs),log10(medianspec(:,currentchannelidx))','r--','LineWidth',1)
            hold on
            plot(log2(freqs),log10(meanspec(:,currentchannelidx))','r','LineWidth',1)
           imagesc(log2(freqs),specpowerbins,(powerdist_best))
            %hold on
            plot(log2(freqs),log10(medianspec(:,currentchannelidx))','r--','LineWidth',1)
            plot(log2(freqs),log10(meanspec(:,currentchannelidx))','r','LineWidth',1)
            LogScale('x',2);LogScale('y',10);
            %caxis([0 0.9])
            axis xy
            axis tight
                    plot(log2(sprange(1)).*[1 1],get(gca,'ylim'),'r--')
            plot(log2(sprange(2)).*[1 1],get(gca,'ylim'),'r--')
         %   colorbar
            xlabel('f (Hz)');ylabel('Power')
            title('Power Density: Best Channel')
            legend('Median Power','Mean Power','location','southwest')
           % ylim([-1 1])

     plotx = linspace(-pi,3*pi,100);      
        subplot(3,2,4)
        hold on
            imagesc(phasebins,powerbins,gammahist)
            imagesc(phasebins+2*pi,powerbins,gammahist)
            plot(gammaangle,powerbins,'.k')
            plot(gammaangle+2*pi,powerbins,'.k')
            plot(plotx,cos(plotx),'k')
            colormap(gca,rwbcolormap2)
            axis xy
            axis tight
            plot([1 1].*angle(spindlegammaPAC),get(gca,'ylim'),'k--')
            plot([1 1].*angle(spindlegammaPAC)+2.*pi,get(gca,'ylim'),'k--')
            ColorbarWithAxis([-0.5 0.5],['Mean High Gamma (',num2str(garange(1)),'-',num2str(garange(2)),') Power (Z(log))'])
            caxis([-0.5 0.5])
            title('Spindle-Gamma Coupling: Best Channel')
          %  xlim([-pi 3*pi]);ylim(powerbins([1 end]))
            xlabel('Spindle Phase');ylabel('Spindle Power (Z(log))')

    subplot(3,3,9)
        gasphist = hist3([spPOWER,gaPOWER],{powerbins,powerbins});
        imagesc(powerbins,powerbins,gasphist)
        axis xy
        xlabel('Spindle Power');ylabel('Gamma Power')

    NiceSave(['SpPower_chan',num2str(currentchannelnumber)],figfolder,filename)


    %% Get pSpindle Intervals
    skewthresh = 0.025;
    skewthresh = 0.03;
    peakthresh = 1.5;
    peakthresh = 1.75;
    durthresh = 0.4;
    
    smoothskew = smooth([0;gammaskew],3);
    smoothskew = smoothskew(2:end);
    
    threshidx = find(smoothskew>=skewthresh,1,'first');
    powerthresh = powerbins(threshidx);
    if isempty(powerthresh); powerthresh=peakthresh; end

    spoverthresh = spPOWER>powerthresh;
    pSpindleInts = IDXtoINT(spoverthresh);
    pSpindleInts = cellfun(@(x) x./sf_best,pSpindleInts,'UniformOutput',false);
    pSpindleInts = pSpindleInts{1};

    pSpDur = pSpindleInts(:,2)-pSpindleInts(:,1);
    pSpInterDur = pSpindleInts(2:end,1)-pSpindleInts(1:end-1,2);

    %%
    viewwin = 0.25;
    pSpLFP_epochs = IsolateEpochs2(spPOWER,pSpindleInts,[viewwin viewwin],sf_best,'includeNaN');

    %% Spindle Features

    numpSp = length(pSpDur);
    nums = length(t_NREMLFP)./sf_best;
    ratepSp = numpSp./nums;
    [~,sortdur] = sort(pSpDur);

    maxPower = cellfun(@(X) max(X([round(viewwin*sf_best):end-round(viewwin*sf_best)])),pSpLFP_epochs);
    [~,sortmaxP] = sort(maxPower);


    %% from pSpindles to Spindles

    SpindleIDX = pSpDur>=durthresh & maxPower>=peakthresh;
    SpindleInts = pSpindleInts(SpindleIDX,:);

    SpDur = SpindleInts(:,2)-SpindleInts(:,1);
    SpInterDur = SpindleInts(2:end,1)-SpindleInts(1:end-1,2);

    numSp = length(SpDur);
    rateSp = numSp./nums;
    [~,sortdur_sp] = sort(SpDur);

    %%
    [pSpNaNMat,t_align_pSp] = CellToNaNMat(pSpLFP_epochs,viewwin,viewwin,sf_best);
    [SpNaNMat,t_align_Sp] = CellToNaNMat(pSpLFP_epochs(SpindleIDX),viewwin,viewwin,sf_best);
    %% Figure: Spindle Thresholds
    sppowermap = makeColorMap([1 1 1],[0 0.5 0],[0.8 0.5 0]);
    figure
    colormap(gcf,sppowermap);
        subplot(4,2,5)
            imagesc(t_align_pSp,[1 length(pSpLFP_epochs)],pSpNaNMat(:,sortdur)')
            xlim([-viewwin 1])
            ColorbarWithAxis([-2 2],'Power (Z(log))')
        subplot(4,2,7)
            imagesc(t_align_Sp,[1 numSp],SpNaNMat(:,sortdur_sp)')
            xlim([-viewwin 1])
            ColorbarWithAxis([-2 2],'Power (Z(log))')
    %     subplot(4,2,3)
    %         imagesc(t_align,[1 length(pSpLFP_epochs)],NaNMat(:,sortmaxP)')
    %         ColorbarWithAxis([-2 2],'Power (Z(log))')
    %         xlim([-viewwin 1])
        subplot(2,2,2)
            plot(pSpDur,maxPower,'.')
            hold on
            plot(get(gca,'xlim'),powerthresh.*[1 1],'k--')
            plot(get(gca,'xlim'),peakthresh.*[1 1],'r--')
            plot(durthresh.*[1 1],get(gca,'ylim'),'r--')
            xlabel('pSpindle Duration (s)');ylabel('Max Power (Z(log))')
            title('From pSpindles to Spindles')


    subplot(4,2,1)
        plot(powerbins,smoothskew,'k','LineWidth',1)
        hold on
        plot([powerbins(1) powerthresh],skewthresh*[1 1],'k--')
        plot(peakthresh.*[1 1],get(gca,'ylim'),'r--')
        plot(powerthresh.*[1 1],[0 skewthresh],'k--') 
        axis tight
        legend('Gamma Coupling','Min Power Thresh (On/Offset)','Max Power Thresh','location','northwest')
        xlabel('Spindle Power (Z(log))');
        ylabel('Gamma PAC (mrl)');
        title('Power Thresholds')

    histbins = linspace(-2,1,30);
    %histbins = linspace(0,2,30);
    pSpDurhist = hist(log10(pSpDur),histbins);
    InterpSpDurhist = hist(log10(pSpInterDur),histbins);
    subplot(4,2,6)
    hold on
        plot(histbins,InterpSpDurhist,'r')
        plot(histbins,pSpDurhist,'k')
        LogScale('x',10)
        axis tight
        xlabel('Duration (s)')
        ylabel('# pSpindles')
        text(-0.25,0.9*max(pSpDurhist),[num2str(round(ratepSp,1)),' pSpindles/s_N_R_E_M'])

    histbins = linspace(-1,2,30);
    SpDurhist = hist(log10(SpDur),histbins);
    InterSpDurhist = hist(log10(SpInterDur),histbins);
    subplot(4,2,8)
    hold on
        plot(histbins,InterSpDurhist,'r')
        plot(histbins,SpDurhist,'k')
        LogScale('x',10)
        axis tight
        xlabel('Duration (s)')
        ylabel('# Spindles')
        text(0.75,0.9*max(SpDurhist),[num2str(round(rateSp,1)),' Spindles/s_N_R_E_M'])

    if numSp~=0
    t_temp = [1:length(NREMLFP_detectchan)]./sf_best;
    randspind = randi(numSp,1);
    randwin =  SpindleInts(randspind,:)+1.5*[-1 1];   
        subplot(4,2,3)

            hold on
            plot(t_temp,spLFP,'r')
            plot(t_temp,spPOWER,'b')
            plot(randwin,powerthresh.*[1 1],'k--')
            plot(randwin,peakthresh.*[1 1],'r--')
            plot(SpindleInts(randspind,:),3.*[1 1],'r','Linewidth',2)
            plot(t_temp,NREMLFP_detectchan,'k','linewidth',1)
            xlim(randwin)
            ylim(4*[-1 1])
    end


    NiceSave(['pSpindleThresh',num2str(currentchannelnumber)],figfolder,filename)

close all  
    %% Mean peak-aligned LFP during Spindles for all Cortical Channels

    %Convert time back into recording time from NREM time
    SpindleIntNREMIDX = round(SpindleInts.*sf_best);
    SpindleInt_realtime{cc} = t_NREMLFP_best(SpindleIntNREMIDX);
    if size(SpindleInt_realtime{cc},2)==1;
        SpindleInt_realtime{cc}=SpindleInt_realtime{cc}';
    end
    Detectedon{cc} = currentchannelnumber.*ones(size(SpindleInt_realtime{cc}(:,1)));

    %% Channel Properties to save
    
    ChannelProperties.powerbins = powerbins;
    ChannelProperties.phasebins = phasebins;
    ChannelProperties.gammaPACmagbypower(:,cc) = gammaskew;
    ChannelProperties.gammaPACanglebypower(:,cc) = gammaangle;
    ChannelProperties.gammapowerhist{cc} = gammahist;
    ChannelProperties.gammaPAC(cc) = spindlegammaPAC;
  
end

%% Join Spindles
minseparation = 0.05;
allspindleints = cat(1,SpindleInt_realtime{:});
alldetectedchannels = cat(1,Detectedon{:});
[mergedspindleints,mergedspindlechannels] = MergeSeparatedInts(allspindleints,minseparation);
mergedspindlechannels = cellfun(@(X) alldetectedchannels(X),mergedspindlechannels,'UniformOutput',false');

numSp = length(mergedspindleints(:,1));
SpDur = mergedspindleints(:,2)-mergedspindleints(:,1);

spindlechannelmat = zeros(numSp,numchannels);
for ss = 1:numSp
    spindlechannelidx = find(ismember(channels2use,mergedspindlechannels{ss}));
    spindlechannelmat(ss,spindlechannelidx) = 1;
end
%%
if numdetectchannels>1

    intoverlap = zeros(numchannels);
    intunderlap = zeros(numchannels);
    intsetoverlap = zeros(numchannels);
    for aa = 1:numchannels
        for bb = aa:numchannels
            [int1overlap,int2overlap]=FindOverlappingInts(SpindleInt_realtime{aa},SpindleInt_realtime{bb});
            intoverlap(aa,bb) = sum(int1overlap)./length(int1overlap);
            intoverlap(bb,aa) = sum(int2overlap)./length(int2overlap);

            intunderlap(aa,bb) = sum(int1overlap==0)./length(int1overlap);
            intunderlap(bb,aa) = sum(int2overlap==0)./length(int2overlap);
            numsp(aa) = length(int1overlap);

            %A and B / A or B
            intsetoverlap(aa,bb) = sum(spindlechannelmat(:,aa)&spindlechannelmat(:,bb))./sum(spindlechannelmat(:,aa)|spindlechannelmat(:,bb));
            intsetoverlap(bb,aa) = intsetoverlap(aa,bb);
        end
    end


    %%
    figure
        plot(mergedspindleints(:,1),(spindlechannelmat),'o')

    imagesc(spindlechannelmat')

    %%
    figure
    hist(sum(spindlechannelmat,2),1:numchannels)
    xlabel('Number of Channels Detected on');ylabel('Number of Spindles')

    %% Figure: Channel Properties
    %add : spindlegammaPAC phase, histogram of #spindles detected on each
    %channel
    gammaPACmagbypower = ChannelProperties.gammaPACmagbypower;
    gammaPACanglebypower = ChannelProperties.gammaPACanglebypower;
    gammacolormap = makeColorMap([1 1 1],[0.8 0 0]);

    figure
    subplot(2,3,2)
    colormap(gca,hsv)
    imagesc(powerbins,[1 numchannels],gammaPACanglebypower')
    alpha((gammaPACmagbypower-skewthresh)'./0.1)
    hold on
    %ColorbarWithAxis(skewthresh+[0 0.1],'Gamma PAC')
    %ColorbarWithAxis([0 0.15],'Gamma PAC')
    title('Spindle-Gamma PAC')
    xlabel('Spindle Power')

    xplot = linspace(-pi,3*pi,50);
    subplot(3,3,8)
        imagesc(xplot,[0 0.02],[xplot xplot])
        colormap(gca,hsv)
        axis xy
        hold on
    plot(angle(ChannelProperties.gammaPAC),abs(ChannelProperties.gammaPAC),'k.')
    plot(angle(ChannelProperties.gammaPAC)+2*pi,abs(ChannelProperties.gammaPAC),'k.')
    plot(xplot,cos(xplot).*0.01+0.01,'k')
    xlim([-pi 3*pi]);ylim([0 max(abs(ChannelProperties.gammaPAC))])
    xlabel('PAC spindle angle');ylabel('PAC mag')


            subplot(2,3,1)
                plot(log2(freqs),bsxfun(@(X,Y) X+Y,(meanskew(:,repchannels)-1).*15,-repchannels),'k','LineWidth',1)
                hold on
                plot(log2(freqs),bsxfun(@(X,Y) X+Y,zeros(size((meanskew(:,repchannels)))),-repchannels+(meanskew(1,repchannels)-1).*15),'--k')
                xlabel('f (Hz)');ylabel('Channel')
                LogScale('x',2)
                axis tight
                title('Spindle Power Skew')
                plot(spindlebandskew(1:end)-1+log2(freqs(end)),-[0:numchannels-1],'ro-')
                plot(spindlebandskew(currentchannelidx)-1+log2(freqs(end)),-[currentchannelidx-1],'r*')
                plot(log2(freqs(end)).*[1 1],get(gca,'ylim'),'k')
                plot(log2(sprange(1)).*[1 1],get(gca,'ylim'),'k--')
                plot(log2(sprange(2)).*[1 1],get(gca,'ylim'),'k--')

        subplot(2,3,3)
            imagesc(intoverlap')
            colorbar
            xlabel('Also detected on channel');
            ylabel('Channel');
            title('Co-Detection')

            subplot(3,3,9)
            colormap(gca,rwbcolormap)
                imagesc(spPowerCorr)
                hold on
                xlabel('Channel');ylabel('Channel')
                set(gca,'Ytick',[]);set(gca,'XTick',[]);
                caxis([-1 1])
                colorbar
                title('Correlation: Spindle Power')

     NiceSave('ChannelProps',figfolder,filename)

end
        
%% Visualization of Spindles - mean and example for each recording
loadwin = 1; %s to load outside of spindles
loadints = bsxfun(@(X,Y) X+Y,mergedspindleints,loadwin.*[-1 1]);
loadints = MergeSeparatedInts(loadints,0);

lfp = [];
indices = [];
for i = 1:length(loadints(:,1)),
    duration = loadints(i,2)-loadints(i,1);
    start = loadints(i,1);
    % Load data
    data = LoadBinary_DL(lfpfilename,'duration',duration,'frequency',Fs,'nchannels',Par.nChannels,'start',start,'channels',(channels2use)+1);
    t = start:(1/Fs):(start+(length(data)-1)/Fs);t=t';
    lfp = [lfp ; t data];
    indices = [indices ; i*ones(size(t))];
end

lfp(:,2:end) = zscore(lfp(:,2:end));

%% Calculate gamma power
[~,gapower_all] = FiltNPhase(lfp(:,2:end),garange,Fs);
gapower_all_norm = bsxfun(@(X,Y) X./Y, gapower_all,mean(gapower_all,1));
gapower_all = zscore(log10(gapower_all));
%% Calculate Spindle
%Filter the reference channel in spindle band and identify peaks
[spLFP_all,spPOWER_all,spPHASE_all] = FiltNPhase(lfp(:,2:end),sprange,Fs);
spPOWER_all = zscore(log10(spPOWER_all));
%% Spindle depth-power, depth-gammaPAC phase
spindlepowerdepth = zeros(numSp,numchannels);
singlespindlePAC = zeros(numSp,numchannels);
for ss = 1:numSp
    [~,inspindleidx] = RestrictInts(lfp(:,1),mergedspindleints(ss,:));
    spindlepowerdepth(ss,:) = mean(spPOWER_all(inspindleidx,:),1);
    [~,maxpeak] = max(spindlepowerdepth(ss,:));
    
    singlespindlePAC(ss,:) = mean(gapower_all_norm(inspindleidx,:).*exp(1i.*spPHASE_all(inspindleidx,:)));
    singlespindlepeakPAC(ss) = singlespindlePAC(ss,maxpeak);
end


%%
[~,sortpower10] = sort(spindlepowerdepth(:,10));

%%
[peakpower,peakchan] = max(spindlepowerdepth,[],2);
%singlespindlepeakPAC = singlespindlePAC(:,peakchan); 
[~,sortpeakchan] = sort(peakchan);
[~,sortpeakpower] = sort(peakpower);
%% Spindles by power depth
% 
% figure
%     subplot(2,2,1)
%     imagesc(spindlepowerdepth(sortpeakchan,:)')
%     ylabel('Channel');xlabel('Spindle (Sorted by peak channel)')
%     colorbar
%     
%     subplot(2,2,2)
%         imagesc(spindlechannelmat(sortpeakchan,:)')
%     ylabel('Channel');xlabel('Spindle (Sorted by peak channel)')
%     colorbar
%     
%     subplot(2,2,3)
%         imagesc(angle(singlespindlePAC(sortpeakchan,:)'))
%         alpha(abs(singlespindlePAC(sortpeakchan,:))'./0.4)
%         colormap(gca,hsv)
%     ylabel('Channel');xlabel('Spindle (Sorted by peak channel)')
%     colorbar
%     
%     subplot(2,2,4)
%         imagesc(abs(singlespindlePAC(sortpeakchan,:)'))
%         
%     ylabel('Channel');xlabel('Spindle (Sorted by peak channel)')
%     colorbar
%     caxis([0 0.2])
%% Mean Spindle/Gamma for spindles detected on each channel
for cc = 1:numdetectchannels
refchannelidx = cc;
refchannelnumber = channels2use(refchannelidx);

onchannelspindles = mergedspindleints(spindlechannelmat(:,refchannelidx)==1,:);
%onchannelspindles = mergedspindleints(peakchan==cc,:);
numonchannelspindles = length(onchannelspindles(:,1));
if numonchannelspindles==0; continue; end

%Identify peaks on the reference channel spindle-band filtered signal
[PKS,LOCS_IDX]= findpeaks(spLFP_all(:,refchannelidx),'MinPeakDistance',0.04.*Fs);
LOCS_t = lfp(LOCS_IDX,1);
%Limit to only peaks in spindles detected on reference channel
LOCS_t = RestrictInts(LOCS_t,onchannelspindles); 


peakwin = 0.2;

%Calculate the mean peak-aligned LFP for all channels with window of 100ms
[peak_epochs] = IsolateEpochs2(lfp(:,2:end),LOCS_IDX,Fs.*[peakwin peakwin],1,'includeNaN');
meanpeakLFP= nanmean(cat(3,peak_epochs{:}),3);

%%
[gapeak_epochs] = IsolateEpochs2(gapower_all,LOCS_IDX,Fs.*[peakwin peakwin],1,'includeNaN');
meanpeakgamma= nanmean(cat(3,gapeak_epochs{:}),3);
%%
t_peaks = linspace(-peakwin,peakwin,length(meanpeakLFP(:,1)));
[~,peakmagsort] = sort(meanpeakLFP(t_peaks==0,:));
%% Mean Spindle Peak

figure
    imagesc(t_peaks([1 end]),[1 numchannels],(meanpeakgamma(:,:))')
    hold on
    colormap(rwbcolormap2)
    plot(t_peaks,bsxfun(@(X,Y) X+Y,-meanpeakLFP(:,repchannels).*5,repchannels),'k','LineWidth',1)
    plot(0,refchannelidx,'ko')
    ColorbarWithAxis(0.1.*[-1 1],'Mean Gamma Power (Z(log))')
    %colorbar
    ylabel('Channel')
    legend('Mean LFP')
    xlabel('t (s) - aligned to reference spindle peak')
    xlim(0.1.*[-1 1])

NiceSave(['MeanPeakAignedSpindle',num2str(refchannelnumber)],figfolder,filename)
%% Example Spindle
viewwin = 0.3;
exspindle = randi(numonchannelspindles,1);
viewint = onchannelspindles(exspindle,:)+viewwin*[-1 1];
t_int = lfp(:,1)>=viewint(1) & lfp(:,1)<=viewint(2);
gammacolormap = makeColorMap([1 1 1],[0.8 0 0]);
spindlecolor = [0 0.6 0];
figure

colormap(gammacolormap)

imagesc(viewint,[1 numchannels],gapower_all(t_int,:)')
hold on
    plot(viewint(1)+viewwin.*[1 1],get(gca,'ylim'),'color',spindlecolor,'LineWidth',2)
    plot(viewint(2)-viewwin.*[1 1],get(gca,'ylim'),'color',spindlecolor,'LineWidth',2)
    plot(lfp(:,1),-spLFP_all(:,refchannelidx)+refchannelidx,'color',spindlecolor,'LineWidth',2)
    plot(lfp(:,1),bsxfun(@(X,Y) X+Y,-lfp(:,repchannels).*0.5,repchannels),'k','LineWidth',1)
    xlim(viewint)
   ColorbarWithAxis([1 3],'Gamma Power (Z(log))')
NiceSave(['ExampleSpindle',num2str(refchannelnumber)],figfolder,filename)

close all
end
%% OUTPUTS

Spindles.ints = mergedspindleints;
%Spindles.timephasemap = 
Spindles.detectedon = Detectedon;
Spindles.detectchannelmat = spindlechannelmat;
%Spindles.gammapowerphase
Spindles.powerdepthprofile = spindlepowerdepth;
Spindles.singlespindlePAC = singlespindlePAC;
Spindles.peakpowerchan = peakchan;
Spindles.peakpowerPAC = singlespindlepeakPAC;


end

