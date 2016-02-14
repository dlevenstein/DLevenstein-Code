function [ structout ] = CollapseStruct( structin,dim,combine )
%structout = CollapseStruct( structin,dim,combine ) Combines elements in a
%structure array
%
%INPUT
%   structin    struct(N).fields structure array with N elements.
%   dim         dimension along which to combine each element in the
%               structure array. (default: 2)
%   (optional)
%       combine     can take 'mean' or 'median' instead of concatenating
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

