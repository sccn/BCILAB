function res = arg_tovals(spec,direct)
% Convert a 'rich' argument report into a 'vals' report.
% Vals = arg_tovals(Rich)
%
% In: 
%   Rich : a 'rich' argument report, as obtained via arg_report('rich',some_function)
%
%   Direct : whether to endow the result with an 'arg_direct' flag set to true, which indicates to 
%            the function taking the Vals struct that the contents of the struct directly correspond
%            to workspace variables of the function. If enabled, contents of Vals must be changed
%            with care - for example, removing/renaming fields will likely lead to errors in the
%            function. (default: true)
%
% Out:
%   Vals : a 'vals' argument report, as obtained via arg_report('vals',some_function) this data
%          structure can be used as a valid argument to some_function.
%
% Examples:
%   % report arguments of myfunction
%   report = arg_report('rich',@myfunction)
%   % convert the report to a valid argument to the function
%   values = arg_tovals(report);
%
% See also:
%   arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-10-18

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if ~exist('direct','var')
    direct = false; end

% remove unassigned specifiers
spec = spec(~strcmp(unassigned,{spec.value}));
% evaluate expressions
expressions = strcmp('expression',{spec.type}) & cellfun('isclass',{spec.value},'char');
if any(expressions)
    try
        [spec(expressions).value] = dealout(evalin('base',format_cellstr({spec(expressions).value}))); 
    catch
        for e=find(expressions)
            try
                spec(e).value = evalin('base',spec(e).value);
            catch
            end
        end
    end
end
% and replace by structs
res = struct('arg_direct',{direct});
for k=1:length(spec)
    if isstruct(spec(k).children)
        % has children: replace by struct
        val = arg_tovals(spec(k).children,direct);
    else
        % no children: take value (and possibly convert to double)
        val = spec(k).value;
        if spec(k).to_double && isinteger(val)
            val = double(val); end
    end
    % and assign the value
    res.(spec(k).names{1}) = val;
end
res.arg_direct = direct;


% like deal(), except that the inputs are given as a cell array instead of a comma-separated list
function varargout = dealout(argin)
varargout = argin;


% format a non-empty cell-string array into a string
function x = format_cellstr(x)
x = ['{' sprintf('%s, ',x{1:end-1}) x{end} '}'];