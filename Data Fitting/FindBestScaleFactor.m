function [ bestsf,dwelltimes_sim_sf,KSSTAT,selectionparm ] = FindBestScaleFactor(dwelltimes_exp,dwelltimes_sim,sf_range)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
% dwelltimes_exp = dwelltimes;
% dwelltimes_sim = dwell;
% sf_range = [0.0001 0.02];
EXFIGS = false;
logtransform = true;
%%
numsf = 300;
if length(sf_range)==1
    scalefactors = sf_range;
    numsf = 1;
elseif length(sf_range)==2
    scalefactors = linspace(sf_range(1),sf_range(2),numsf);
else
    error('wrong # scale factors')
end
%%
if any(isnan(dwelltimes_sim.UP)) || any(isnan(dwelltimes_sim.DOWN))
    KSSTAT.DOWN = nan; KSSTAT.UP = nan; selectionparm = nan;
    dwelltimes_sim_sf.UP = nan;  dwelltimes_sim_sf.DOWN = nan; 
    bestsf = nan;
    return
end

if logtransform
    dwelltimes_exp.UP = log10(dwelltimes_exp.UP);
    dwelltimes_exp.DOWN = log10(dwelltimes_exp.DOWN);
end

%Loop through possible scaling factors and do the KSTest
for ss = 1:numsf

    sf_test = scalefactors(ss);
    dwelltimes_sim_test.DOWN = dwelltimes_sim.DOWN.*sf_test;
    dwelltimes_sim_test.UP = dwelltimes_sim.UP.*sf_test;
    
    if logtransform
        dwelltimes_sim_test.UP = log10(dwelltimes_sim_test.UP);
        dwelltimes_sim_test.DOWN = log10(dwelltimes_sim_test.DOWN);
    end
    
    [~,~,KSSTAT.DOWN(ss)] = kstest2(dwelltimes_sim_test.DOWN,dwelltimes_exp.DOWN);
    [~,~,KSSTAT.UP(ss)] = kstest2(dwelltimes_sim_test.UP,dwelltimes_exp.UP);
    
end

selectionparm = (1-KSSTAT.UP).*(1-KSSTAT.DOWN);
[~,bestsfidx] = max(selectionparm);

bestsf = scalefactors(bestsfidx);
dwelltimes_sim_sf.DOWN = dwelltimes_sim.DOWN.*bestsf;
dwelltimes_sim_sf.UP = dwelltimes_sim.UP.*bestsf;


%% Panels for Figures and all that

%Bad sf's...
badsf1.DOWN = dwelltimes_sim_sf.DOWN.*2;
badsf1.UP = dwelltimes_sim_sf.UP.*2;
badsf2.DOWN = dwelltimes_sim_sf.DOWN./2;
badsf2.UP = dwelltimes_sim_sf.UP./2;

    
if logtransform
    dwelltimes_sim_sf_plot.UP = log10(dwelltimes_sim_sf.UP);
    dwelltimes_sim_sf_plot.DOWN = log10(dwelltimes_sim_sf.DOWN);
    
    badsf1.DOWN = log10(badsf1.DOWN);
    badsf1.UP = log10(badsf1.UP);
    badsf2.DOWN = log10(badsf2.DOWN);
    badsf2.UP = log10(badsf2.UP);
end


if EXFIGS
figure
    subplot(2,1,1)
    plot(scalefactors,KSSTAT.DOWN,'r','LineWidth',1)
    hold on
    plot(scalefactors,KSSTAT.UP,'g','LineWidth',1)
    plot(scalefactors,selectionparm,'k','LineWidth',2)
    legend('K-S: DOWN','K-S: UP','Similarity','location','southeast')
    xlabel('Time Scale Factor (s/tau_r)')
    ylabel('K-S Statistic / Similarity Metric')
    title('Finding Time Scale Factor')
    subplot(2,1,2)
hold on
plot(sort(dwelltimes_sim_sf_plot.UP),(1:length(dwelltimes_sim_sf_plot.UP))./length(dwelltimes_sim_sf_plot.UP),'g','LineWidth',2)
plot(sort(dwelltimes_sim_sf_plot.DOWN),(1:length(dwelltimes_sim_sf_plot.DOWN))./length(dwelltimes_sim_sf_plot.DOWN),'r','LineWidth',2)
plot(sort(dwelltimes_exp.UP),(1:length(dwelltimes_exp.UP))./length(dwelltimes_exp.UP),'g--','LineWidth',2)
plot(sort(dwelltimes_exp.DOWN),(1:length(dwelltimes_exp.DOWN))./length(dwelltimes_exp.DOWN),'r--','LineWidth',2)
legend('UP: Simulation','DOWN: Simulation','UP: Experiment','DOWN: Experiment',...
    'location','southeast')
title(['Best Scale Factor: ',num2str(round(bestsf,4))])
xlabel('Dwell Time (s)');ylabel('Cumulative')


%% FOR FIGURES: Similarity Metric Map
ks_DOWN = linspace(0,1,100);
ks_UP = linspace(0,1,100);
simmetric = (1-ks_DOWN)'*(1-ks_UP);

figure
imagesc(ks_UP,ks_DOWN,simmetric)
colorbar
xlabel('K-S DOWN');ylabel('K-S UP')
ColorbarWithAxis([0 1],'Similarity')
axis xy

%% FOR POSTER: Example
ff = figure;
    subplot(8,4,1)
    histogram(dwelltimes_sim_sf_plot.UP,15,'facecolor','none')
    hold on
    histogram(dwelltimes_sim_sf_plot.DOWN,15,'facecolor','g','FaceAlpha',0.2)
    axis tight
    box off
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    xlabel('In Silico Dwell Times (tau_r)')
        subplot(8,4,5)
    histogram(dwelltimes_exp.UP,20,'facecolor','none')
    hold on
    histogram(dwelltimes_exp.DOWN,20,'facecolor','g','FaceAlpha',0.2)
    axis tight
    xlim([-1.5 1.5])
    box off
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    xlabel('In Vivo Dwell Times (s)')
    
    subplot(2,6,3:4)
    plot(scalefactors,KSSTAT.DOWN,'g','LineWidth',1)
    hold on
    plot(scalefactors,KSSTAT.UP,'k','LineWidth',1)
    plot(scalefactors,selectionparm,'color',[0.9 0.6 0],'LineWidth',2)
    plot([bestsf bestsf],[0 selectionparm(bestsfidx)],'k','LineWidth',0.1)
    plot([0 bestsf],[selectionparm(bestsfidx) selectionparm(bestsfidx)],'k','LineWidth',0.1)
    legend('K-S: DOWN','K-S: UP','Similarity','location','southoutside')
    xlabel('Time Scale Factor (s/tau_r)')
    ylabel({'K-S Distance', '/ Similarity'})
   % title('Finding Time Scale Factor')
    set(gca,'ytick',[0 1]);
    xlim([0 sf_range(end)])
    set(gca,'xtick',get(gca,'xlim'));
    
    
fitcolormap = makeColorMap([1 1 1],[0 0.6 0],[0.9 0.6 0]) ;
subplot(5,4,9)
%axes('position',[0.55 0.8 0.1 0.1])
imagesc(ks_UP,ks_DOWN,simmetric)
colorbar
colormap(fitcolormap)
xlabel('K-S_D_O_W_N');ylabel('K-S_U_P')
ColorbarWithAxis([0.1 0.9],'Similarity')
axis xy
set(gca,'xtick',[0 1]);set(gca,'ytick',[0 1]);    
    
    subplot(4,3,3)
    hold on
plot(sort(dwelltimes_sim_sf_plot.UP),(1:length(dwelltimes_sim_sf_plot.UP))./length(dwelltimes_sim_sf_plot.UP),'k--','LineWidth',2)
plot(sort(dwelltimes_sim_sf_plot.DOWN),(1:length(dwelltimes_sim_sf_plot.DOWN))./length(dwelltimes_sim_sf_plot.DOWN),'g--','LineWidth',2)
plot(sort(dwelltimes_exp.UP),(1:length(dwelltimes_exp.UP))./length(dwelltimes_exp.UP),'k','LineWidth',2)
plot(sort(dwelltimes_exp.DOWN),(1:length(dwelltimes_exp.DOWN))./length(dwelltimes_exp.DOWN),'g','LineWidth',2)
%title(['Best Scale Factor: ',num2str(round(bestsf,4))])
xlabel('Dwell Time (s)');ylabel('Cumulative')
box on
axis tight
xlim([-1.5 1.5])
set(gca,'ytick',[0 1]);
%set(gca,'xtick',[0 5 10]);

    subplot(4,3,6)
    hold on
plot(sort(dwelltimes_sim_sf_plot.UP),(1:length(dwelltimes_sim_sf_plot.UP))./length(dwelltimes_sim_sf_plot.UP),'k--','LineWidth',2)
plot(sort(dwelltimes_sim_sf_plot.DOWN),(1:length(dwelltimes_sim_sf_plot.DOWN))./length(dwelltimes_sim_sf_plot.DOWN),'g--','LineWidth',2)
plot(sort(dwelltimes_exp.UP),(1:length(dwelltimes_exp.UP))./length(dwelltimes_exp.UP),'k','LineWidth',2)
plot(sort(dwelltimes_exp.DOWN),(1:length(dwelltimes_exp.DOWN))./length(dwelltimes_exp.DOWN),'g','LineWidth',2)
legend('UP: Simulation','DOWN: Simulation','UP: Experiment','DOWN: Experiment',...
    'location','southeast')
%title(['Best Scale Factor: ',num2str(round(bestsf,4))])
xlabel('Dwell Time (s)');ylabel('Cumulative')
box on
axis tight
xlim([-1.5 1.5])
set(gca,'ytick',[0 1]);
%set(gca,'xtick',[0 5 10]);


%CDFs for bad time scale factors
    subplot(4,3,9)
    hold on
plot(sort(badsf1.UP),(1:length(dwelltimes_sim_sf_plot.UP))./length(dwelltimes_sim_sf_plot.UP),'k--','LineWidth',2)
plot(sort(badsf1.DOWN),(1:length(dwelltimes_sim_sf_plot.DOWN))./length(dwelltimes_sim_sf_plot.DOWN),'g--','LineWidth',2)
plot(sort(dwelltimes_exp.UP),(1:length(dwelltimes_exp.UP))./length(dwelltimes_exp.UP),'k','LineWidth',2)
plot(sort(dwelltimes_exp.DOWN),(1:length(dwelltimes_exp.DOWN))./length(dwelltimes_exp.DOWN),'g','LineWidth',2)
%legend('UP: Simulation','DOWN: Simulation','UP: Experiment','DOWN: Experiment',...
%    'location','southeast')
%title(['Best Scale Factor: ',num2str(round(bestsf,4))])
xlabel('Dwell Time (s)');ylabel('Cumulative')
box on
axis tight
xlim([-1.5 1.5])
set(gca,'ytick',[0 1]);
%set(gca,'xtick',[0 5 10]);

    subplot(4,3,12)
    hold on
plot(sort(badsf2.UP),(1:length(dwelltimes_sim_sf_plot.UP))./length(dwelltimes_sim_sf_plot.UP),'k--','LineWidth',2)
plot(sort(badsf2.DOWN),(1:length(dwelltimes_sim_sf_plot.DOWN))./length(dwelltimes_sim_sf_plot.DOWN),'g--','LineWidth',2)
plot(sort(dwelltimes_exp.UP),(1:length(dwelltimes_exp.UP))./length(dwelltimes_exp.UP),'k','LineWidth',2)
plot(sort(dwelltimes_exp.DOWN),(1:length(dwelltimes_exp.DOWN))./length(dwelltimes_exp.DOWN),'g','LineWidth',2)
%legend('UP: Simulation','DOWN: Simulation','UP: Experiment','DOWN: Experiment',...
%    'location','southeast')
%title(['Best Scale Factor: ',num2str(round(bestsf,4))])
xlabel('Dwell Time (s)');ylabel('Cumulative')
box on
axis tight
xlim([-1.5 1.5])
set(gca,'ytick',[0 1]);
%set(gca,'xtick',[0 5 10]);


% 
% 
% 
set(ff,'PaperOrientation','landscape');
set(ff,'PaperUnits','normalized');
set(ff,'PaperPosition', [0 0 1 1]);
saveas(ff,'FitFigure','pdf')


end
    
%% Outputs
KSSTAT.DOWN = KSSTAT.DOWN(bestsfidx);
KSSTAT.UP = KSSTAT.UP(bestsfidx);
selectionparm = selectionparm(bestsfidx);
end


