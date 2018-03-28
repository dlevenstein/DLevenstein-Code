function [ R,parms,bins,phaseamprate ] = ParametricLFPCoupling( spiketimes,LFP,frange,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%INPUT
%   spiketimes          {Ncells} cell array of spiketimes for each cell
%                                   or TSObject
%   LFP
%   (optional)
%       'sf_LFP'
%       'waveparms'...
%       'nfreqs' %not implemented yet......
%       'int'
%       'downsample'
%
%COUPLING PARMS
%   alpha: power modulation
%   phi: preferred phase
%   thmod: phase-coupling magnitude
%   pthresh: power threshold for phase coupling
%
%
%DLevenstein Feb 2017
%% inputParse for Optional Inputs and Defaults
p = inputParser;

defaultSf = 1250;

defaultInt = [0 Inf];
checkInt = @(x) size(x,2)==2 && isnumeric(x) || isa(x,'intervalSet');

defaultNfreqs = 1;
defaultNcyc = 5;
checkFrange = @(x) isnumeric(x) && length(x(1,:)) == 2 && length(x(:,1)) == 1;

defaultDOWN = false;

addParameter(p,'sf_LFP',defaultSf,@isnumeric)
addParameter(p,'int',defaultInt,checkInt)
addParameter(p,'nfreqs',defaultNfreqs,@isnumeric)
addParameter(p,'ncyc',defaultNcyc,@isnumeric)
addParameter(p,'DOWNSAMPLE',defaultDOWN,@isnumeric)

parse(p,varargin{:})
%Clean up this junk...
sf_LFP = p.Results.sf_LFP;
ints = p.Results.int;
nfreqs = p.Results.nfreqs;
ncyc = p.Results.ncyc;

%% Deal with input types

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
    for pp = 1:numpop
        if isempty(spiketimes{pp})
            spiketimes{pp} = {};
            continue
        end
        for cc = 1:length(spiketimes{pp})
            spiketimestemp{cc} = Range(spiketimes{pp}{cc},'s');
        end
        spiketimes{pp} = spiketimestemp;
        clear spiketimestemp
    end
    spiketimes = cat(2,spiketimes{:});
    numcells = length(spiketimes);
elseif isa(spiketimes,'cell')
    numcells = length(spiketimes);
end

if isa(ints,'intervalSet')
    ints = [Start(ints,'s'), End(ints,'s')];
end

t_LFP = (1:length(LFP))'/sf_LFP;


%% Get only spike times in ints


for cc = 1:numcells
    [spiketimes{cc},~] = RestrictInts(spiketimes{cc},ints);
    spiketimes{cc} = spiketimes{cc}(:,1);
end

[~,LFP_intidx] = RestrictInts(t_LFP,ints);

%% Filter the LFP

%Filter
[~,amp,phase] = FiltNPhase(LFP,frange,sf_LFP);

%Normalize Power
amp = NormToInt(amp,ints,sf_LFP,'percentile');


%%

%Get Power/Phase at each spike
spikeamp = cellfun(@(X) interp1(t_LFP,amp,X,'nearest'),spiketimes,'uniformoutput',false);
spikephase = cellfun(@(X) interp1(t_LFP,phase,X,'nearest'),spiketimes,'uniformoutput',false);

%Get Power/Phase Time Occupancy
t_amp = interp1(t_LFP,amp,t_LFP(LFP_intidx),'nearest');
t_phase = interp1(t_LFP,phase,t_LFP(LFP_intidx),'nearest');

%%


numbins = 20;

%Power/Phase Bins
phasebins = linspace(-pi,pi,numbins+1);phasebins=phasebins(1:end-1)+0.5.*diff(phasebins(1:2));
powerbins = linspace(0,1,numbins+1);powerbins=powerbins(1:end-1)+0.5.*diff(powerbins(1:2));


%Calculate Power/Phase Spike Histogram
phaseamphist = cellfun(@(X,Y) hist3([X,Y],{powerbins,phasebins}),spikeamp,spikephase,'UniformOutput',false);
%Calculate Power/Phase Time Histogram
phaseamphist_t = hist3([t_amp,t_phase],{powerbins,phasebins});
phaseamphist_t = phaseamphist_t./sf_LFP;

%Normalize Rate
totaltime = sum(diff(ints,1,2));
meanrate = cellfun(@(X) length(X)./totaltime,spiketimes,'UniformOutput',false);
phaseamprate = cellfun(@(X) X./phaseamphist_t,phaseamphist,'UniformOutput',false);
phaseamprate = cellfun(@(X,Y) X./Y,phaseamprate,meanrate,'UniformOutput',false);

%Rho/Theta for input to fitting function
rho = repmat(powerbins',1,numbins);
theta = repmat(phasebins,numbins,1);


%%

ft = fittype( 'ParametricPowerPhase( rho,theta,alpha,beta,phi,pthresh,thmod )','independent',{'rho','theta'} );

options = fitoptions(ft);
options.StartPoint = [0,0,0,0.5,0];
options.Lower = [-Inf 0 -pi 0 0];
options.Upper = [Inf Inf 3*pi 1 Inf];
options.MaxIter = 1e15;
options.MaxFunEvals = 1e15;

couplingfit = cellfun(@(X) fit( [rho(:),theta(:)], X(:), ft, options),phaseamprate,'UniformOutput',false);

%%
couplingparms = cellfun(@(X) coeffvalues(X),couplingfit,'UniformOutput',false);
parmsconf = cellfun(@(X) confint(X),couplingfit,'UniformOutput',false);
%%
for cc = 1:numcells
    R{cc} = ParametricPowerPhase( rho,theta,...
    couplingparms{cc}(1),couplingparms{cc}(2),couplingparms{cc}(3),couplingparms{cc}(4),couplingparms{cc}(5));
end

couplingparms = cat(1,couplingparms{:});
parms.alpha = couplingparms(:,1);
parms.beta = couplingparms(:,2);
parms.phi = couplingparms(:,3);
parms.pthresh = couplingparms(:,4);
parms.thmod = couplingparms(:,5);

bins.phasebins = phasebins;
bins.powerbins = powerbins;

end

