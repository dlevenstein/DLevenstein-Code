function [int1overlap,int2overlap,int2whichint1,int1whichint2] = FindOverlappingInts(ints1,ints2)
%[int1overlap,int2overlap] = FindOverlappingInts(ints1,ints2) finds
%intervals in sets ints1 and ints2 that overlap in time. Returns logicals
%that indicate which intervals in the sets have an interval in the other
%set that is at the same time.
%%
if isempty(ints1) && isempty(ints2)
    int1overlap = [];
    int2overlap = [];
    return
elseif isempty(ints1)
    int1overlap = [];
    int2overlap = zeros(size(ints2(:,1)));
    return
elseif isempty(ints2)
    int1overlap = zeros(size(ints1(:,1)));
    int2overlap = [];
    return
end

numints1 = length(ints1(:,1));
numints2 = length(ints2(:,1));

int1overlap = false(numints1,1);
int2overlap = false(numints2,1);
for ii = 1:numints1
    beforestart = ints2<=ints1(ii,1);
    aroundstart = beforestart(:,1)~=beforestart(:,2);
    int2overlap = int2overlap|aroundstart;
    afterend = ints2>=ints1(ii,2);
    aroundend = afterend(:,1)~=afterend(:,2);
    int2overlap = int2overlap|aroundend;
    betweenstartend = all([~beforestart,~afterend],2);
    int2overlap = int2overlap|betweenstartend;
    
    %which intervals are the ones that overlap this one?
    whichoverlap = [find(aroundend) find(aroundstart) find(betweenstartend)];
    whichoverlap = unique(whichoverlap);
    if isempty(whichoverlap)
        whichoverlap = [];
    end
    int2whichint1{ii} = whichoverlap;
end

for ii = 1:numints2
    beforestart = ints1<=ints2(ii,1);
    aroundstart = beforestart(:,1)~=beforestart(:,2);
    int1overlap = int1overlap|aroundstart;
    afterend = ints1>=ints2(ii,2);
    aroundend = afterend(:,1)~=afterend(:,2);
    int1overlap = int1overlap|aroundend;
    betweenstartend = all([~beforestart,~afterend],2);
    int1overlap = int1overlap|betweenstartend;
    
    %which intervals are the ones that overlap this one?
    whichoverlap = [find(aroundend) find(aroundstart) find(betweenstartend)];
    whichoverlap = unique(whichoverlap);
    if isempty(whichoverlap)
        whichoverlap = [];
    end
    int1whichint2{ii} = whichoverlap;
end
1;

