function [dwelltimes_sim_sf,dwellfit] = WCadapt_DwellTimeMatch(dwelltimestomatch,simdwelltimes,scalefactor)
%[matcheddwelltimes] = WCadapt_DwellTimeMatch(dwelltimestomatch,simdwelltimes)
%
%INPUT
%   dwelltimestomatch   structure of UP/DOWN dwell times
%       .UP
%       .DOWN
%   simdwelltimes
%   scalefactor (optional) single value or range
%   
%%
if isempty(dwelltimestomatch.UP) || isempty(dwelltimestomatch.DOWN)
    dwelltimes_sim_sf = [];dwellfit=[];
    return
end

if ~exist('scalefactor','var')
    scalefactor = [0.001 0.025];
end

numsims = length(simdwelltimes);
parfor_progress(numsims);
tstart = tic;
parfor nn = 1:numsims

    timespent=toc(tstart);
    percdone = parfor_progress;

    estimatedtotal = timespent./(percdone./100);
    estimatedremaining = estimatedtotal-timespent;
    display(['Percent Done: ',num2str(percdone),...
        '.  Time Spent: ',num2str(round(timespent./60,1)),...
        '.  Est. Total Time: ',num2str(round(estimatedtotal./60,1)),...
        'min.  ETR: ',num2str(round(estimatedremaining./60,1)),'min.'])
    
        [dwellfit(nn).bestsf,dwelltimes_sim_sf(nn),...
            dwellfit(nn).KSSTAT,dwellfit(nn).selectionparm]...
            = FindBestScaleFactor(dwelltimestomatch,simdwelltimes(nn),scalefactor);
end

parfor_progress(0);

end

