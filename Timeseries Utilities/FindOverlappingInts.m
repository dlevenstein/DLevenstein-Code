function [ int1overlap, int2overlap ] = FindOverlappingInts( ints1,ints2 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here6
%% Test
% testints1 = [0 3; 5 6; 7 10; 13 15; 14.4 14.5];
% testints2 = [1 4; 4.5 7; 12 13; 14 15];
% ints1 = testints1;
% ints2 = testints2;

%%
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
end

