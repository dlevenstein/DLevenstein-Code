function [] = LogScale( whichaxis,logbase)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if strcmp(whichaxis,'y') || strcmp(whichaxis,'xy')
    range = get(gca,'YLim');
    range(1) = floor(range(1));range(2) = ceil(range(2));
    set(gca,'YTick',[range(1):range(2)])
    set(gca,'YTickLabels',logbase.^[range(1):range(2)])
end

if strcmp(whichaxis,'x') || strcmp(whichaxis,'xy')
    range = get(gca,'XLim');
    range(1) = floor(range(1));range(2) = ceil(range(2));
    set(gca,'XTick',[range(1):range(2)])
    set(gca,'XTickLabels',logbase.^[range(1):range(2)])
end

end

