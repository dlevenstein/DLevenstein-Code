function [ spikesynch ] = SpikeSynch( spiketimes,synchwin,synchtype )
%[spikesynch] = SpikeSynch(spiketimes,synchwin) calculates population 
%synchrony of all other cells within a synchwin window around each spike in
%spiketimes
%
%INPUT
%   spiketimes  {ncells} cell array of [nspikes] spiketimes
%   synchwin    size of window around each spike to count (s)
%   synchtype   'numcells' or 'numspikes' (default numcells)
%               'numcells' will be in units of proportion of pop (i.e. 0-1)
%               'numspikes' will be in units of spks/cell/s (Hz)
%OUPUTS
%   spikesynch  population synchrony in synchwin window around
%               each spike in spiketimes
%
%DLevenstein Winter 2017
%%
numcells = length(spiketimes);

if ~exist('synchtype','var') 
    synchtype = 'numcells';
end

for cc = 1:numcells
    spikesynch{cc} = zeros(size(spiketimes{cc}));
%cc = 1;
display(['Cell: ',num2str(cc),' of ',num2str(numcells)])


%%MAKE THIS INTO A FUNCTION AND USE CELLFUN??%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numspks = length(spiketimes{cc});
    othercells = setdiff(1:numcells,cc);
    for kk = 1:numspks
        %ss = 2;
        spkwin = spiketimes{cc}(kk) + 0.5*synchwin*[-1 1];
        %All othercell spikes within window
        inwin = cellfun(@(X) X>=spkwin(1) & X<=spkwin(2),spiketimes(othercells),'UniformOutput',false);
        switch synchtype
            case 'numcells'
                synchinwin = sum(cellfun(@(X) max(X),inwin));
            case 'numspikes'
                synchinwin = sum(cellfun(@(X) sum(X),inwin))./synchwin;
        end
        spikesynch{cc}(kk) = (synchinwin)./(numcells-1);
    end  
%%MAKE THIS INTO A FUNCTION AND USE CELLFUN??%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Alternatively
%join all spikes into a [spiketime cellnum] vector, might be faster

end


end

