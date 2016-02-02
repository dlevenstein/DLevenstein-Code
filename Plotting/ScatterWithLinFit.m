function [  ] = ScatterWithLinFit( x,y,color,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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


[R,p] = corr(x,y);%,'type','Spearman');



%%
props.color = color;
plot(x,y,'.','Color',props.color,'markersize',10)
lsline

xbouns = get(gca,'Xlim');
ybouns = get(gca,'ylim');
textloc = [xbouns(1),ybouns(2) - 0.1*(ybouns(2)-ybouns(1))];
text(textloc(1),textloc(2),{['R = ',num2str(R)],['p = ',num2str(p)]},'color',props.color)
xlim(xbouns);ylim(ybouns);
end

