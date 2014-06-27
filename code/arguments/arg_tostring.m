function str = arg_tostring(val,strip_direct,indent,subindent)
% Convert an argument value to a string.
% String = arg_tostring(Value,StripDirect,Indent,SubIndent)
%
% Converts a value to a string that is functionally equivalent to the result of hlp_tostring, but
% which follows the formatting guidelines of the arg functions (particularly arg_subtoggle and
% arg_subswitch), plus support indentation and line breaks for readability.
%
% In:
%   Value : the value to convert; this is a result as produced by arg_report for
%           either the 'rich', 'lean', or 'nvps' report type.
%
%   StripDirect : whether to strip occurrences of arg_direct from the value
%                 (default: false)
%
%   Indent : initial indentation of the string (default: 0)
%
%   SubIndent : indentation increment for nested lines (default: 4)
%
% Out:
%   String : a human-readable, pretty-printed string that should yield a value under eval()
%            that is functionally equivalent as the original Value (as far as the argument system
%            is concerned).
%
% See also:
%   arg_report, hlp_tostring
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-25

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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
if nargin < 2
    strip_direct = false; end
if nargin < 3
    indent = 0; end
if nargin < 4
    subindent = 4; end

% if this is an argument specification, turn it into a cell array of 
% name-value pairs with human-readable names
if all(isfield(val,{'head','children','alternatives'}))
    val = arg_tovals(val,[],'HumanReadableCell'); end

str = '';           % the output string
selection = [];     % value of the arg_selection entry, if any

if iscellstr(val(1:2:end)) && mod(length(val),2) == 0
    % got name-value pairs: special formatting
    if strip_direct
        % optionally get rid of any arg_direct flags
        p = find(strcmp('arg_direct',val(1:2:end)));
        if ~isempty(p)
            val([2*p-1 2*p]) = []; end
    end
    p = find(strcmp('arg_selection',val(1:2:end)));
    if ~isempty(p)
        % got an argument with an arg_selection element: special print formatting applies
        % first we strip the arg_selection element (generally)
        selection = val{2*p(end)};
        val([2*p-1 2*p]) = [];
        % implement formatting shortcuts
        if isequal(selection,false)
            % we have a disabled subtoggle argument: replaced by 'off'
            str = hlp_tostring('off');
        elseif isempty(val)
            if isequal(selection,true)
                % enabled arg_subtoggle with no overridden values: replaced by 'on'
                str = hlp_tostring('on');
            elseif ischar(selection)
                % arg_subswitch with no overridden values: replaced by the selection string
                str = hlp_tostring(selection);
            end
        end
    end
    % if the string has not yet been assigned, generate it now
    if isempty(str)
        shift = @(k) repmat(' ',1,k);
        if ischar(selection)
            str = ['{' '''' selection '''' ' ...'];
        elseif isequal(selection,true) || isequal(selection,[])
            str = '{ ...';
        else
            error(['Unsupported arg_selection form: ' hlp_tostring(selection)]);
        end
        % pretty-print name-value pairs with line breaks
        for c=1:2:length(val)
            str = [str sprintf('\n') shift(indent+subindent) '''' val{c} '''' ...
                ', ' arg_tostring(val{c+1},strip_direct,indent+subindent,subindent) ' ...'];  %#ok<AGROW>
        end 
        str = [str(1:end-4) '}'];
    end
elseif islogical(val) && isscalar(val)
    % these two values are not printed in a particularly pretty way by hlp_tostring
    if val
        str = 'true';
    else
        str = 'false';
    end
else
    % any other value (opaque as far as the arg system is concerned)
    str = hlp_tostring(val);
end
