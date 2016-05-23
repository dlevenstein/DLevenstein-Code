function [epochcell] = TimeNormalize(epochcell,tbins)
%[epochcell] = TimeNormalize(epochcell,tbins) time normalizes all epochs in 
% epochcell to be the same number of time bins. If epochs are 2d (example:
% spectrogram) must have time on first dimension.
%
%INPUTS
%   -epochcell  [Nepochs x 1] cell array of time epochs, 
%               each is [Nt x Nvar]
%   -tbins      number of bins (length in time) the output epochs should be
%
%OUTPUTs
%   -epochcell  time normalized
%
%Last Updated: 10/8/15
%DLevenstein
%%


epochlengths = cellfun(@(A) length(A(:,1)),epochcell);
numepochs = length(epochcell);


newepoch = [1 1];    %Placeholder
for ee = 1:numepochs;
    display(['Epoch: ',num2str(ee),' of ',num2str(numepochs)])
    resamplefact = tbins/epochlengths(ee);
    tol = 0.0001;
    while length(newepoch(:,1)) ~= tbins
        [P,Q] = rat(resamplefact,tol);
        if P==0
            tol = tol/10;
            continue
        end
        if P*Q >=2^20 | tol<1e-300  %Avoid crashing resample...
            epochcell{ee}([1,end],:) = [];
            epochlengths(ee) = epochlengths(ee)-2;
            resamplefact = tbins/epochlengths(ee);
            tol = 0.0001;
            continue
        end
        newepoch = resample(epochcell{ee},P,Q);
        tol = tol/10;
    end
    epochcell{ee} = newepoch;
    newepoch = [1 1];
end


end

