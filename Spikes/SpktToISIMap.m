function [ ISImap,ISIdensity,ISIhist,histbins ] = SpktToISIMap(spiketimes,varargin )
%[ ISImap,ISIdensity,ISIhist,autocorr ] = SpktToISIMap( spiketimes,int )
%returns the ISImap and other ISI statistics from the spike times of a set
%of cells, with the option to only look at spikes in a given interval
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] vectors containing 
%               spiketimes for each cell
%   int         [Nints x 2] set of intervals within which to consider
%               spikes. (need to make this optional)
%
%OUTPUT
%   
%
%
%DLevenstein 2016
%%
%spiketimes = Se;
%int = [End(SWSPacketInts,'s')-20 End(SWSPacketInts,'s')];
%int = REMInts;
%int = StateIntervals.Spindles;
%% Optional input arguments
p = inputParser;

defaultInts = [0 Inf];
checkInts = @(x) size(x,2) == 2 || isa(x,'intervalSet');

defaultSHOWFIG = false;

defaultNumbins = 30;


addParameter(p,'showfig',defaultSHOWFIG,@islogical)
addParameter(p,'ints',defaultInts,checkInts)
addParameter(p,'numbins',defaultNumbins,@isnumeric)

parse(p,varargin{:})
%Clean up this junk...
SHOWFIG = p.Results.showfig; 
int = p.Results.ints;
numbins = p.Results.numbins;
%%
numcells = length(spiketimes);
%% Deal with Input Types
if numcells == 0
    display('case no cells needs outputs')
    pause
    return
end

%Spiketimes can be: tsdArray of cells, cell array of cells, cell array of
%tsdArrays (multiple populations)
if isa(spiketimes,'tsdArray')
    numcells = length(spiketimes);
    for cc = 1:numcells
        spiketimestemp{cc} = Range(spiketimes{cc},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
elseif isa(spiketimes,'cell') && isa(spiketimes{1},'tsdArray')
    numpop = length(spiketimes);
    lastpopnum = 0;
    for pp = 1:numpop
        if length(spiketimes{pp})==0
            spiketimes{pp} = {};
            popcellind{pp} = [];
            continue
        end
        for cc = 1:length(spiketimes{pp})
            spiketimestemp{cc} = Range(spiketimes{pp}{cc},'s');
        end
        spiketimes{pp} = spiketimestemp;
        popcellind{pp} = [1:length(spiketimes{pp})]+lastpopnum;
        lastpopnum = popcellind{pp}(end);
        clear spiketimestemp
    end
    spiketimes = cat(2,spiketimes{:});
    numcells = length(spiketimes);
    subpop = 'done';
end

if isa(int,'intervalSet')
    int = [Start(int,'s'), End(int,'s')];
end


%%
ISIs = cellfun(@diff,spiketimes,'UniformOutput',false);
ISIt = cellfun(@(X) X(2:end-1),spiketimes,'UniformOutput',false);
ISIn = cellfun(@(X) X(1:end-1),ISIs,'UniformOutput',false);
ISInp1 = cellfun(@(X) X(2:end),ISIs,'UniformOutput',false);

%% Restrict to Ints
[ISIt,keepidx] = cellfun(@(X) RestrictInts(X,int),ISIt,'UniformOutput',false);
ISIn = cellfun(@(X,Y) X(Y),ISIn,keepidx,'UniformOutput',false);
ISInp1 = cellfun(@(X,Y) X(Y),ISInp1,keepidx,'UniformOutput',false);
ISImap = cellfun(@(X,Y) [X Y],ISIn,ISInp1,'UniformOutput',false);
%%
histbins = linspace(-2.5,1.5,numbins);
ISIhist = cellfun(@(X) hist(log10(X),histbins),ISIn,'UniformOutput',false);
ISIdensity = cellfun(@(X,Y) hist3(log10([X,Y]),{histbins,histbins}),ISIn,ISInp1,'UniformOutput',false);
ISIdensity = cat(3,ISIdensity{:});
ISIhist = cat(1,ISIhist{:});

%%
if SHOWFIG
    numspikes = cellfun(@length,ISIn);
    [~,sortrate] = sort(numspikes);
    
    for cc = 1:numcells
    cellnum = sortrate(cc);
    subidx = mod(cc-1,10)+1;
    histoffset = 20;
    if subidx == 1
    figure
    end
    if subidx > 5
        subidx = subidx+10;
        histoffset = 10;
    end
       subplot(5,5,subidx)
            imagesc(histbins,histbins,ISIdensity{cellnum})
            axis xy
            hold on
            plot(log10(1/10.*[1 1]),get(gca,'ylim'),'r--')
            plot(log10(1/20.*[1 1]),get(gca,'ylim'),'r--')
            plot(get(gca,'xlim'),log10(1/10.*[1 1]),'r--')
            plot(get(gca,'xlim'),log10(1/20.*[1 1]),'r--')
            LogScale('xy',10)
            ylabel('ISI n+1 (s)')
            xlim(histbins([1 end])+[-0.05 0.05]);ylim(histbins([1 end])+[-0.05 0.05])
       subplot(5,5,5+subidx)
            plot(log10(ISIn{cellnum}),log10(ISInp1{cellnum}),'.')
            hold on
            plot(log10(1/10.*[1 1]),get(gca,'ylim'),'r--')
            plot(log10(1/20.*[1 1]),get(gca,'ylim'),'r--')
            plot(get(gca,'xlim'),log10(1/10.*[1 1]),'r--')
            plot(get(gca,'xlim'),log10(1/20.*[1 1]),'r--')
            LogScale('xy',10)
            ylabel('ISI n+1 (s)')
            xlim(histbins([1 end])+[-0.05 0.05]);ylim(histbins([1 end])+[-0.05 0.05])
       subplot(10,5,histoffset+subidx)
            bar(histbins,ISIhist{cellnum})
            hold on
            plot(log10(1/10.*[1 1]),get(gca,'ylim'),'r--')
            plot(log10(1/20.*[1 1]),get(gca,'ylim'),'r--')
            plot(get(gca,'xlim'),log10(1/10.*[1 1]),'r--')
            plot(get(gca,'xlim'),log10(1/20.*[1 1]),'r--')
            LogScale('x',10)
            xlim(histbins([1 end])+[-0.05 0.05]);
            xlabel('ISI n (s)');ylabel('ISI n+1 (s)')
            ylim([0 max(ISIhist{cellnum})])

    end      
end

end

