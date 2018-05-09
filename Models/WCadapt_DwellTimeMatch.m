function [dwelltimes_sim_sf,dwellfit] = WCadapt_DwellTimeMatch(dwelltimestomatch,simdwelltimes,varargin)
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
p = inputParser;
addParameter(p,'scalefactor',[0.001 0.025])
addParameter(p,'showprogress',true,@islogical)

parse(p,varargin{:})
scalefactor = p.Results.scalefactor;
showprogress = p.Results.showprogress;

%%
if isempty(dwelltimestomatch.UP) || isempty(dwelltimestomatch.DOWN)
    dwelltimes_sim_sf = [];dwellfit=[];
    return
end



numsims = length(simdwelltimes);
if showprogress
        parfor_progress(numsims);
end
    tstart = tic;


parfor nn = 1:numsims

    if showprogress
        timespent=toc(tstart);
        percdone = parfor_progress;

        estimatedtotal = timespent./(percdone./100);
        estimatedremaining = estimatedtotal-timespent;
        display(['Percent Done: ',num2str(percdone),...
            '.  Time Spent: ',num2str(round(timespent./60,1)),...
            '.  Est. Total Time: ',num2str(round(estimatedtotal./60,1)),...
            'min.  ETR: ',num2str(round(estimatedremaining./60,1)),'min.'])
    end
    
        [dwellfit(nn).bestsf,dwelltimes_sim_sf(nn),...
            dwellfit(nn).KSSTAT,dwellfit(nn).selectionparm]...
            = FindBestScaleFactor(dwelltimestomatch,simdwelltimes(nn),scalefactor);
end

if showprogress
    parfor_progress(0);
end

end

