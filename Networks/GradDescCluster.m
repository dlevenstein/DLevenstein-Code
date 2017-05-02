function [ cluass ] = GradDescCluster( cohmat )
%SpatialGradDescCluster()
%
%INPUT
%   
%
%% DEV
%load('/Users/dlevenstein/Dropbox/Research/Buzsaki Lab/AmygdalaHPC/LFPSiteMapAnalysis/gacohState.mat')
%cohmat = corrmat.sws;
%% Parms and defaults
numsteps = 500000;

numinit = 20; %number of initial random clusters


%% Initiate and Run
numsites = size(cohmat,1);
%Start from random initial assignments
cluass = randi(numinit,numsites,1);

for ss = 1:numsteps
    cluass = DescOneStep(cluass,cohmat);
end


%% FUNCTION: OneStep
    function [newcluass] = DescOneStep(oldcluass,cohmat)
        clus = unique(oldcluass);
        numclus = length(clus);

        %Pick a random site and calculate the energy for it's current cluster
        sitepick = randi(length(oldcluass),1);
        %othersites = setdiff(1:length(oldcluass),sitepick);
        clu_A = oldcluass(sitepick);
        N_A = sum(oldcluass == clu_A);
        E_A = -(1./N_A).*sum(cohmat(sitepick,oldcluass==clu_A));

        %Calculate energy gap for switching site to other clusters
        energygap = zeros(size(clus));
        for cc = 1:numclus
            clu_B = clus(cc);
            if clu_B==clu_A
                continue
            end
            N_B = sum(oldcluass == clu_B)+1; 
            %note: site i included in A and B
            E_B = -(1./N_B).*(sum(cohmat(sitepick,[find(oldcluass==clu_B);sitepick])));
            energygap(cc) = E_A - E_B;
        end

        %Switch site to cluster with largest energy gap
        newcluass = oldcluass;
        newcluass(sitepick) = clus(energygap==max(energygap));

    end
end

