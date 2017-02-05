function [correlation,signif] = ScatterWithLinFit( x,y,color,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Options/Defaults

sig = 'p0';
corrtype = 'spearman';


%%
if length(x(1,:))>length(x(:,1))
    x = x';
end
if length(y(1,:))>length(y(:,1))
    y = y';
end

nonans = ~isnan(x) & ~isnan(y);
x = x(nonans);
y = y(nonans);

noinfs = ~isinf(x) & ~isinf(y);
x = x(noinfs);
y = y(noinfs);

if sum(ismember(varargin, 'RemoveOutlierX'))
    [~,xoutlierind] = max(abs(x-mean(x)));
    x(xoutlierind) = [];
    y(xoutlierind) = [];
end

%Correlation
switch corrtype
    case 'pearson'
        [R,p] = corr(x,y);
        corrtext = ['R = ',num2str(round(R,3))];
        correlation = R; signif = p;
    case 'spearman'
        [rho,p] = corr(x,y,'type','Spearman');
         corrtext = ['\rho = ',num2str(round(rho,2))];
         correlation = rho; signif = p;
end



switch sig
    case 'p0'
        sigtext = ['   p = ',num2str(round(p,3))];
        
    case 'conf'
        %Confidence Interval
        [pf, S] = polyfit(x,y,1);
        % yfit = polyval(pf,ex);
        % NB: ... "ex" is a logspace bunch of values fitting the range of the data
        estY=polyval(pf,x); % the estimate of the 'y' value for each point.
        SlopeCI = polyparci(pf,S,0.95);

        corrtext = ['R = ',num2str(pf(1))];
        sigtext = [' +/-',num2str(diff(SlopeCI(:,1))/2)];
        
    case 'none'
        sigtext = [];
end

%%
props.color = color;
plot(x,y,'.','Color',props.color,'markersize',10)
axis tight
lsline

xbouns = get(gca,'Xlim');
ybouns = get(gca,'ylim');
textloc = [xbouns(1),ybouns(2) - 0.1*(ybouns(2)-ybouns(1))];
textloc = [xbouns(1)+ 0.05*(xbouns(2)-xbouns(1)),ybouns(2) + 0.05*(ybouns(2)-ybouns(1))];
text(textloc(1),textloc(2),{[corrtext,sigtext]},'color',props.color)
xlim(xbouns);ylim(ybouns);
end

