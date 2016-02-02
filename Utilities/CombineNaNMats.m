function [ bigmat ] = CombineNaNMats( mats )
%UNTITLED2 Summary of this function goes here
%   mats: cell array of differently sized trial matrices from multiple
%   sessions or something

bigmat = [];
for r = 1:length(mats)
    onerec = mats{r};
   % onerec = (onerec-nanmean(onerec(:)))./nanstd(onerec(:));
    
    lengthpack = length(onerec(1,:));
    if r==1
        bigmat = onerec;
        continue
    elseif lengthpack > length(bigmat(1,:))
        lengthall = length(bigmat(1,:));
        widthall = length(bigmat(:,1));
        newnans = lengthpack - lengthall;
        bigmat = [bigmat NaN(widthall,newnans)];
    else
        lengthall = length(bigmat(1,:));
        widthall = length(onerec(:,1));
        newnans = lengthall - lengthpack;  
        onerec = [onerec NaN(widthall,newnans)];
    end
    bigmat = vertcat(bigmat,onerec);
end


end

