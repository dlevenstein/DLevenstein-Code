function [ crange ] = SpecColorRange( spec )
%[ crange ] = SpecColorRange( spec ) calculates a good color range for a
%spectrogram and sets the plot to the good range.
%
%DLevenstein 2017
%%

[~,mu,sig] = zscore(spec);
crange = [min(mu)-2.5*max(sig) max(mu)+2.5*max(sig)];
caxis(crange)

end

