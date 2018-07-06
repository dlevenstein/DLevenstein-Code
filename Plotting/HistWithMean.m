function [histmean,histstd,histp] = HistWithMean(data,varargin)
%[histmean,histstd,histp] = HistWithMean(data,numbins)
%   Detailed explanation goes here
%%
p = inputParser;
addParameter(p,'numbins',10);
addParameter(p,'color','k');
addParameter(p,'showtext',true);


parse(p,varargin{:})

numbins = p.Results.numbins;
color4face = p.Results.color;
showtext = p.Results.showtext;



noinfs = ~isinf(data);
data = data(noinfs);

if ~exist('color4face','var')
    color4face = [1 1 1];
end


histmean = nanmean(data);
histstd = nanstd(data);


[~,histp]=ttest(data);

hold on
histogram(data,numbins,'facecolor',color4face)

histvals = hist(data,numbins);
ylims = get(gca,'ylim');
meany = ylims(2)*1.01;
meany = max(histvals)*1.25;


plot(histmean,meany,'r+')
plot(histmean+[-histstd,histstd],meany*[1,1],'k-','LineWidth',1)

if showtext
xbouns = get(gca,'Xlim');
ybouns = get(gca,'ylim');
textloc = [xbouns(1),ybouns(2) - 0.05*(ybouns(2)-ybouns(1))];
text(textloc(1),textloc(2),{['mean = ',num2str(histmean)],['p = ',num2str(histp)]})
xlim(xbouns);ylim(ybouns);
end
end

