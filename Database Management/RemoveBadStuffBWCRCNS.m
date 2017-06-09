function [  ] = RemoveBadStuffBWCRCNS( basepath )

basepath = pwd;
basename = bz_BasenameFromBasepath(basepath);

mkdir(fullfile(basepath,'ClusBeforeOverwrite'))

%% Blank out spikes from bad shanks 
load(fullfile(basepath,[basename '_BasicMetaData.mat']))
totalnumberofshanks = length(bmd.Par.SpkGrps);

badshanks = setdiff(1:totalnumberofshanks,bmd.goodshanks);
for bidx = 1:length(badshanks)
    tshank = badshanks(bidx);
    cluname = fullfile(basepath,[basename '.clu.' num2str(tshank)]);
    
    clu = load(cluname);
    clu = zeros(size(clu));%blank out all to noise cluster
    clu(1) = 1;%add denotation that there is 1 cluster here (first entry to .clu is 1 now)
    
    movefile(cluname,fullfile(basepath,'ClusBeforeOverwrite',[basename '.clu.' num2str(tshank)]));

    %write to new clu
    fid = fopen(cluname,'w'); 
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid
end

%% from good shanks, find badcells, call them noise
load(fullfile(basepath,[basename '_SStable.mat']),'badcells')
load(fullfile(basepath,[basename '_SAll.mat']),'shank','cellIx')

shankstoload = unique(shank(badcells.allbadcells));

for sidx = 1:length(shankstoload)
    tshank = shankstoload(sidx);
    cluname = fullfile(basepath,[basename '.clu.' num2str(tshank)]);
    clu = load(cluname);
    clu = clu(2:end); % toss the first sample to match res/spk files
    
    movefile(cluname,fullfile(basepath,'ClusBeforeOverwrite',[basename '.clu.' num2str(tshank)]));
    
    badsinshank = (badcells.allbadcells([tshank == shank(badcells.allbadcells)]));%overall number of this cell (ie number considering all shanks)
    for bidx = length(badsinshank):-1:1
        thiscellix = cellIx(badsinshank(bidx));%idx in this shank
        clu(clu==thiscellix) = 0;
        %scoot down the number of all clusters above the one deleted
        clu(clu> thiscellix) = clu(clu> thiscellix)-1;
    end
    
    clu = cat(1,length(unique(clu)),clu);
    %rewrite this clu
    fid=fopen(cluname,'w'); 
    fprintf(fid,'%.0f\n',clu);
    fclose(fid);
    clear fid
end

end

