function [ spiketimes ] = ShuffleSpikelabels(spiketimes)
%[ spiketimes ] = ShuffleSpikelabels(spiketimes) jitters spiketimes
%of all cells within a designated window.
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] vectors of spiketimes for
%               each cell
%               -or-
%               [Nspikes x 2] array of [spiketimes cellnumber]
%
%OUTPUT
%   spiketimes  jittered spiketimes
%
%
%DLevenstein 2017
%%
if iscell(spiketimes) || isa(spiketimes,'tsdArray')
    inputtype = 'cellofcells';
elseif isnumeric(spiketimes) && size(spiketimes,2)==2
    inputtype = 'timecellarray';
else
    error('Your spiketimes input has some issues, eh?')
end

switch inputtype
    case 'cellofcells'
        
        numcells = length(spiketimes);
        if isa(spiketimes,'tsdArray')  %If spikes are in a TSObject
            for cc = 1:numcells
                spiketimestemp{cc} = Range(spiketimes{cc},'s');
            end
            spiketimes = spiketimestemp;
            clear spiketimestemp
        end

        %Count how many spikes each cell should get at the end
        numcellspks = cellfun(@length,spiketimes);

        %Put all the spike times together into a single vector, then shuffle them
        allspks = cat(1,spiketimes{:});
        allspks = allspks(randperm(length(allspks)));

        %Put the shuffled times back into a cell array with same number of spikes
        %per cell
        spiketimes = mat2cell(allspks,numcellspks)';
        spiketimes = cellfun(@sort,spiketimes,'UniformOutput',false);
        
        
    case 'timecellarray'
        
        %Randomly permute cell ID labels
        numspks = length(spiketimes(:,1));
        spiketimes(:,2) = spiketimes(randperm(numspks),2);
        
        
end


end

