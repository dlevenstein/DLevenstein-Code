function [dwelltimes_sim_sf,dwellfit] = WCadapt_DwellTimeMatch(dwelltimestomatch,simdwelltimes)
%[matcheddwelltimes] = WCadapt_DwellTimeMatch(dwelltimestomatch,simdwelltimes)
%
%INPUT
%   dwelltimestomatch   structure of UP/DOWN dwell times
%       .UP
%       .DOWN
%   simdwelltimes
%   
%%


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
            = FindBestScaleFactor(dwelltimestomatch,simdwelltimes(nn),[0.001 0.025]);
end

parfor_progress(0);

end

