function [ structout ] = CollapseStruct( structin,dim,combine,NEST )
%structout = CollapseStruct( structin,dim,combine ) Combines elements in a
%structure array
%
%INPUT
%   structin    struct(N).fields structure array with N elements where each
%               of the N elements has the same fields and field structure.
%               structin can have nested cell arrays or stuctures, but be
%               careful with this... may not work for >1-dimensional
%               structures/cell arrays
%   dim         dimension along which to combine each element in the
%               structure array. (default: 2)
%   (optional)
%       combine     can take 'mean' or 'median' instead of concatenating
%                   or 'std'
%
%       
%% A place for input options
if ~exist('NEST','var')
    NEST = false;
end

if ~exist('combine','var')
    combine = 'justcat';
end

if ~exist('dim','var')
    dim = 2;
end

%%
fields = fieldnames(structin);

%%

for ff = 1:length(fields)
    currentfield = fields{ff};
   % structout = setfield(structout,currentfield,value)
%     if exist('combine','var')
%         if strcmp(combine,'mean')
%             for ii = 1:length(structin)
%                 structin(ii).(currentfield) = nanmean(structin(ii).(currentfield),dim);
%             end
%         elseif strcmp(combine,'median')
%             for ii = 1:length(structin)
%                 structin(ii).(currentfield) = nanmedian(structin(ii).(currentfield),dim);
%             end
%         end
%     end

	if isstruct(structin(1).(currentfield)) & NEST %For Nested Structures
       structout.(currentfield) = cat(1,structin(:).(currentfield));
       structout.(currentfield) = CollapseStruct(structout.(currentfield),dim,combine,true);
       continue
    elseif iscell(structin(1).(currentfield)) & NEST %For cell array in field
        structout.(currentfield) = cat(1,structin(:).(currentfield));
        structout.(currentfield) = cat(dim,structout.(currentfield){:});
    elseif (isstring(structin(1).(currentfield))||ischar(structin(1).(currentfield))) & NEST %For string in field
        structout.(currentfield) = {structin(:).(currentfield)};
    else %For simple array in field
        try
        structout.(currentfield) = cat(dim,structin(:).(currentfield));
        catch
            keyboard
            continue
        end
    end

    
    
    switch combine
        case 'mean'
            structout.(currentfield) = nanmean(structout.(currentfield),dim);
        case 'median'
            structout.(currentfield) = nanmedian(structout.(currentfield),dim);
        case 'std'
            structout.(currentfield) = nanstd(structout.(currentfield),[],dim);
        case 'sem'
            structout.(currentfield) =  nanstd(structout.(currentfield),[],dim)./sqrt(sum(~isnan(structout.(currentfield)),dim));
    end
end



end

