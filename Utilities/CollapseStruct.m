function [ structout ] = CollapseStruct( structin,dim,combine )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
%
%%
fields = fieldnames(structin);

%%
if ~exist('dim','var')
    dim = 2;
end

for ff = 1:length(fields)
    currentfield = fields{ff};
   % structout = setfield(structout,currentfield,value)
    if exist('combine','var')
        if strcmp(combine,'mean')
            for ii = 1:length(structin)
                structin(ii).(currentfield) = nanmean(structin(ii).(currentfield),dim);
            end
        elseif strcmp(combine,'median')
            for ii = 1:length(structin)
                structin(ii).(currentfield) = nanmedian(structin(ii).(currentfield),dim);
            end
        end
    end
    structout.(currentfield) = cat(dim,structin(:).(currentfield));
end



end

