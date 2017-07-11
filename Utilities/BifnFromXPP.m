function [ bifnlines ] = BifnFromXPP( filename )
%[ bifnlines ] = BifnFromXPP( filename )
%%

Ibifn = importdata(filename);

biflineidx = unique(Ibifn(:,4));
bifline2idx = unique(Ibifn(:,5));

bifnlines = Ibifn(:,1);
for bb = biflineidx'
    for bb2 = bifline2idx'
        bifnlines = [bifnlines Ibifn(:,3)];
        bifnlines(Ibifn(:,4)~=bb | Ibifn(:,5)~=bb2,end) =nan;
        
    end
        

end

