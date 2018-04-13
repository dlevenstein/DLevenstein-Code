function [ bifnlines ] = BifnFromXPP( filename )
%[ bifnlines ] = BifnFromXPP( filename )
%%
if ~exist(filename,'file')
    error('filename given does not exist!')
end


Ibifn = importdata(filename);

% biflineidx = unique(Ibifn(:,4));
% bifline2idx = unique(Ibifn(:,5));
% 
% bifnlines = Ibifn(:,1);
% for bb = biflineidx'
%     for bb2 = bifline2idx'
%         bifnlines = [bifnlines Ibifn(:,3)];
%         
%         bifnlines(~(Ibifn(:,4)==bb & Ibifn(:,5)==bb2),end) =nan;
%         
%     end
%         
% 
% end
%%
thisline = 1;
thiscol = 2;
thispos = 1;
thisentry = Ibifn(thispos,[4 5]);
idx(1,1)=1;
for rr = 1:length(Ibifn)
    %If a new line
    if ~isequal(Ibifn(rr,[4 5]),thisentry)
        idx(thisline,2) = rr-1;
        idx(thisline+1,1) = rr;
        thisline= thisline+1;
        thisentry = Ibifn(rr,[4 5]);
    end
end
idx(thisline,2) = rr;

numhopf = 0;
bifnlines = nan(length(Ibifn(:,1)),thisline+1);
bifnlines(:,1) = Ibifn(:,1);
for ll = 1:thisline
    bifnlines(idx(ll,1):idx(ll,2),ll+1) = Ibifn(idx(ll,1):idx(ll,2),3);
    
    %For limit cycles
    if ~isequal(Ibifn(idx(ll,1):idx(ll,2),2),Ibifn(idx(ll,1):idx(ll,2),3))
        numhopf = numhopf+1;
        bifnlines(:,thisline+numhopf) = nan;
        bifnlines(idx(ll,1):idx(ll,2),thisline+numhopf) = Ibifn(idx(ll,1):idx(ll,2),2);
    end
end
    

%% Detect/Remove large jumps
jumpsize = diff(bifnlines);  
jumpthresh = 0.05;  %should set based on jump sizes.....
%figure
%hist(log10(abs(diffsize(:,3))));
%LogScale('x',10)
bifnlines(abs(jumpsize)>jumpthresh)=nan;
%%
%figure
%plot(Ibifn(:,1),Ibifn(:,2),'k.')
