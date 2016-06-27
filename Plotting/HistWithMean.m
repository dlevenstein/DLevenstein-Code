function [histmean,histstd,histp] = HistWithMean(data,numbins)
%[histmean,histstd,histp] = HistWithMean(data,numbins)
%   Detailed explanation goes here
%%

noinfs = ~isinf(data);
data = data(noinfs);



histmean = nanmean(data);
histstd = nanstd(data);


[~,histp]=ttest(data);

hold on
hist(data,numbins)

histvals = hist(data,numbins);
ylims = get(gca,'ylim');
meany = ylims(2)*1.05;
meany = max(histvals)*1.1;


plot(histmean,meany,'r+')
plot(histmean+[-histstd,histstd],meany*[1,1],'k-','LineWidth',1)

xbouns = get(gca,'Xlim');
ybouns = get(gca,'ylim');
textloc = [xbouns(1),ybouns(2) - 0.1*(ybouns(2)-ybouns(1))];
text(textloc(1),textloc(2),{['mean = ',num2str(histmean)],['p = ',num2str(histp)]})
xlim(xbouns);ylim(ybouns);
end

