function res = arg_tovals(spec,set_direct,format,mandatory_check,unassigned_check,expression_check,conversion_check)
% Reformat an argument specification for a function into a valid input argument.
% Vals = arg_tovals(Spec)
%
% In: 
%   Specification : an argument specification (struct array), e.g., from arg_report('lean',...) or
%                   arg_report('rich',...)
%
%   SetDirect : whether to endow the result with an 'arg_direct' flag set to true, which indicates to 
%               the function taking the Vals struct that the contents of the struct directly correspond
%               to workspace variables of the function. If enabled, contents of Vals must be changed
%               with care - for example, removing/renaming fields will likely lead to errors in the
%               function. If set to [], no change will be made. (default: false)
%
%   Format : Output format, can be one of the following:
%             * 'struct' : generate a struct with values assigned to fields (default)
%             * 'cell' : generate a cell array of name-value pairs, using code names 
%             * 'HumanReadableCell' : generate a cell array with human-readable names
%
%   MandatoryCheck : whether to check and produce errors for unspecified mandatory arguments
%                    (default: false)
%
%   UnassignedCheck : whether to check and strip unassigned and skippable arguments from the result
%                     (default: true)
%
%   ExpressionCheck : whether to check and evalute expression arguments from the result
%                     (default: true)
%
%   ConversionCheck : whether to check for type conversions and apply them to the result
%                     (default: true)
%
% Out:
%   Vals : a data structure (struct or cell array of name-value pairs) that can be passed to the
%          function for which the Specification was generated as an input.
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

    % check inputs
    if nargin < 7
        conversion_check = true; 
        if nargin < 6
            expression_check = true; end
        if nargin < 5
            unassigned_check = true; end
        if nargin < 4
            mandatory_check = false; end
        if nargin < 3
            format = 'struct'; end
        if nargin < 2
            set_direct = false; end
    end
    
    % remove skippable arguments
    spec([spec.skippable]) = [];
    
    % get the values
    values = {spec.value};

    % generate errors for mandatory arguments that were not assigned
    if mandatory_check
        mandatory_entries = find(strcmp('__arg_mandatory__',values));
        if mandatory_entries
            error(['The arguments ' format_cellstr({spec(mandatory_entries).first_name}) ' were unspecified but are mandatory.']); end
    end

    % evaluate any expression-typed arguments
    if expression_check
        expressions = find(strcmp('expression',{spec.type}));
        if ~isempty(expressions)
            expressions(~cellfun('isclass',values(expressions),'char')) = [];
            if ~isempty(expressions)
                try
                    [values{expressions}] = celldeal(evalin('base',format_cellstr(values(expressions))));
                catch %#ok<CTCH>
                    for e=expressions
                        try
                            values{e} = evalin('base',values{e});
                        catch %#ok<CTCH>
                        end
                    end
                end
            end
        end
    end

    % remove unassigned arguments
    if unassigned_check
        toprune = strcmp('__arg_unassigned__',values);
        if any(toprune)
            spec(toprune) = [];
            values(toprune) = [];
        end
    end

    % convert values to double if necessary
    if conversion_check
        for k=find([spec.to_double])
            if isinteger(values{k})
                values{k} = double(values{k}); end
        end
    end

    % recursively process structs
    for k=find(~cellfun('isempty',{spec.children}))
        values{k} = arg_tovals(spec(k).children,set_direct,format,mandatory_check,unassigned_check,expression_check,conversion_check); end

    % build the output data structure
    switch format
        case 'struct'
            res = cell2struct(values,{spec.first_name},2);
            if ~isempty(set_direct)
                try
                    res.arg_direct = set_direct; 
                catch %#ok<CTCH>
                    res = struct('arg_direct',set_direct);
                end
            end
        case 'cell'
            res = reshape([{spec.first_name};values],1,[]);
            if ~isempty(set_direct)
                res = [res {'arg_direct',set_direct}]; end
        case 'HumanReadableCell'
            res = reshape([cellfun(@(n)n{min(2,length(n))},{spec.names},'UniformOutput',false);values],1,[]);
            if ~isempty(set_direct)
                res = [res {'arg_direct',set_direct}]; end
        otherwise
            error(['Unrecognized output format: ' hlp_tostring(format)]);
    end            
end


% format a non-empty cell-string array into a string
function x = format_cellstr(x)
    x = ['{' sprintf('%s, ',x{1:end-1}) x{end} '}'];
end
