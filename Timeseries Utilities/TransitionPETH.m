function [epochs,PETH,droppedints] = TransitionPETH(data,int1,int2,window,sf,specparms,tol)
%[epochs,PETH] = TransitionPETH(data,int1,int2,window,sf,specparms)
%
%INPUTS
%   data    (nt x ndim)
%   int1    preceeding interval set
%   int2    postceeding interval set
%          (optional) int1 can be a cell array of interval sets and int2 can
%           be a cell array of names, will compute TransitionPETH for each
%           combination of intervals
%           (optional) int1/2 can be 'all'
%   window  [secondsbefore,secondsafter] time window around transitions
%   sf      sampling frequency of data
%   specparms (optional) if want to calculate wavelet spectrum of data
%       .frange
%       .nfreqs
%       .ncyc   number of cycles for wavelet
%       .down   downsample factor after wavelets
%   
%
%
%TO DO
%    -default specparms so don't have to include them every time
%   -use inputparser - check matlab doc
%
%
%Last Updated: 11/18/15
%DLevenstein
%%
if ~exist('tol','var')
    tol = 1;
end

if isempty(data)
    display('No Data')
    epochs = {};
    PETH = [];
    return
end



%Self-calling loop for transitions between multiple states
if iscell(int1) && iscell(int2)
    numintsets = length(int1);
    for ii = 1:numintsets
        for jj = 1:numintsets
            transname = strcat(int2{ii},int2{jj});
            [epochs.(transname),PETH.(transname)] = TransitionPETH(data,...
                int1{ii},int1{jj},window,sf,specparms);
        end
    end
    return
end

%In case input is in form of TSObjects
if isa(int1,'intervalSet')
    int1 = [Start(int1,'s'), End(int1,'s')];
end
if isa(int2,'intervalSet')
    int2 = [Start(int2,'s'), End(int2,'s')];
end


%Find Transitions from int1 to int2
[transitions1,transitions2] = FindIntsNextToInts(int1,int2,tol);
if isempty(transitions1) && isempty(transitions2)
    epochs = {};
    PETH = [];
    display('No Interval Transitions')
    return
elseif isempty(transitions1)
    transitions = transitions2(:,1);
else
    transitions = transitions1(:,2);
end


%Get the data around state transitions
[epochs,droppedints] = IsolateEpochs2(data,transitions,window,sf);

%For Specta
if exist('specparms','var')
    frange = specparms.frange;
    nfreqs = specparms.nfreqs;
    ncyc = specparms.ncyc;
    space = 'log';
    [freqs,~,epochs] = WaveSpec(epochs,frange,nfreqs,ncyc,1/sf,space);
    epochs = cellfun(@(X) (abs(X')),epochs,'UniformOutput',false);
    
    PETH.freqs = freqs;
    specdown = specparms.down;
    epochs = cellfun(@(X) downsample(X,specdown),epochs,'UniformOutput',false);
end
%%
meandim = ndims(epochs{1})+1;
PETH.mean = nanmean(cat(meandim,epochs{:}),meandim);
% PETH.std = std(cat(meandim,epochs{:}),[],meandim);
% PETH.fano = PETH.std./PETH.mean;
PETH.t = [-window(1):1/sf:window(2)];


%%
if exist('specparms','var')
    figure
        subplot(2,1,1)
            hold on
            imagesc(PETH.t,log2(freqs),PETH.mean')
            plot([0 0],get(gca,'ylim'),'k')
            xlim(PETH.t([1,end]));ylim(log2(freqs([1,end])))
            colorbar
            axis xy
            LogScale('y',2)
            xlabel('t - Aligned to Transition')
            ylabel('f (Hz)')
else
    figure
        hold on
        plot(PETH.t,PETH.mean)
        plot([0 0],get(gca,'ylim'),'k')
        xlabel('t - Aligned to Transition')
end

end

