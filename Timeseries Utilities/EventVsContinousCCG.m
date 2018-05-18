function [CCGmean,t_lag,CCGstd,CCGevents,skippedevents ] = EventVsContinousCCG( data,data_t,events,lagwin )
%Data should be continuous time series associated with time stamps data_t
%events should be in seconds, as should lagwin
%
%INPUTS
%   data
%%
%Find frames closest to WHon, NWHon transitions
event_ints = interp1(data_t,data_t,events,'nearest');
%Get average pupil, EMG aligned to transitions (20s lag)

%estimate dt, number of samples around events to grab
data_dt = mean(diff(data_t));
%t_lag = 20; %s
t_lag_dt = round(lagwin./data_dt);
t_lag = [-t_lag_dt:t_lag_dt]'.*data_dt;

if isempty(events)
    CCGmean = zeros(size(t_lag));
    return
end

%Get the time points closest to events, and samples in window around them
[events_dt] = find(ismember(data_t,event_ints));
events_dt = bsxfun(@(A,B) A+B,events_dt,t_lag_dt.*[-1 1]);

skippedevents = (isnan(event_ints));

[CCGevents,~] = IsolateEpochs2(data,events_dt,0,1,'includeNaN');
CCGevents = cat(2,CCGevents{:});
CCGmean = nanmean(CCGevents,2);
CCGstd = nanstd(CCGevents,[],2);

end

