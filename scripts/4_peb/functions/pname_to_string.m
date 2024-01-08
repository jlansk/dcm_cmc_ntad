function [out,parts] = pname_to_string(pname, region_names, input_names)
% Translates a PEB parameter string e.g. B(1,2,3) to a friendly descriptor
% e.g. B(From region 1 to region 3 input 3)
%
% pname        - parameter name string from PEB e.g. A{2}(1,2,3)
% region_names - cell array of region names
% input_names  - cell array of input (condition) names
%
% out          - friendly name for the parameter
% parts        - cell array for row, col and (fMRI) task input

str = ['(?<field>[A-Za-z0-9\{\},]+)\('... % Match field and open bracket
       '(?<row>\d+)(,|\))'...     % Match row and open bracket or comma
       '(?<col>\d+)?(,|\))?'...   % Match column and open bracket or comma
       '(?<input>\d+)?(,|\))?'];  % Match input and open bracket or comma

parts = regexp(pname, str, 'names');

parts.row   = str2double(parts.row);
parts.col   = str2double(parts.col);
parts.input = str2double(parts.input);

if isempty(region_names)
    out = pname;
    return;
end

% Skip G-matrix for CMC model
if strcmp(parts.field,'G')
    out = pname;
    return;
end

out = [parts.field '-matrix '];
if isnan(parts.col)
    % Row only
    out = sprintf('%s %s', out, region_names{parts.row});
elseif strcmp(parts.field,'C')
    % Row and region
    out = sprintf('%s %s',region_names{parts.row}, input_names{parts.col});
else
    % Row and col
    out = sprintf('%s from %s to %s', ...
        out, ...
        region_names{parts.col}, ...
        region_names{parts.row});
end

if ~isnan(parts.input)
    out = sprintf('%s (%s)', out, input_names{parts.input});
end