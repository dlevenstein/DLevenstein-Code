function [x,y,u,v] = MeanVelocityFromPosition(xpos,ypos,int,sf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
%
%TO DO:
%   -Include confidence.... variability in velocity etc?
%%



if exist('sf','var')
    h = 1/sf;
else
    h = 1;
end


% %Numerical Differentiation using symmetric difference quotient method.
% dx = (xpos(3:end)-xpos(1:end-2))/(2*h);
% dy = (ypos(3:end)-ypos(1:end-2))/(2*h);
% xpos = xpos(2:end-1);
% ypos = ypos(2:end-1);
% tlen = length(xpos);

%Numerical Differentiation using Newton's Method (first order)
dx = diff(xpos)./h;
dy = diff(ypos)./h;
xpos = xpos(1:end-1);
ypos = ypos(1:end-1);
tlen = length(xpos);


if exist('int','var')
    if exist('sf','var')
        int = int*sf;
        int(:,1) = floor(int(:,1))+1;
        int(:,2) = ceil(int(:,2))+1;
    end
    intIDX = INTtoIDX({int},tlen);
    intIDX = logical(intIDX);
    %Restrict X,dX to intervals
    xpos = xpos(intIDX);
    ypos = ypos(intIDX);
    dx = dx(intIDX);
    dy = dy(intIDX);
end
    


x = [];
y = [];
u = [];
v = [];
xvals = unique(xpos);
yvals = unique(ypos);
for xind = 1:length(xvals)
    xx = xvals(xind);
    display(['x: ',num2str(xx),' of ',num2str(xvals(end))])
    for yind = 1:length(yvals)
        yy = yvals(yind);
        x = [x xx];
        y = [y yy];
        
        meandx = mean(dx(xpos==xx & ypos==yy));
        meandy = mean(dy(xpos==xx & ypos==yy));
        u = [u meandx];
        v = [v meandy];
    end
end

end

