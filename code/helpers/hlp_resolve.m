function outargs = hlp_resolve(x,default,context)
% Look up a symbol from the current hlp_scope scope.
% Data = hlp_resolve(Symbol,Default,Context)
%
% More precisely, perform a lookup for a symbol in the current (dynamic) scope formed by (possibly
% nested) hlp_scope(s). Returns a cell array of all the values associated with the symbol. Used to
% implement dynamic scoping across functions, used by parts of the expression & online system.
%
% In:
%   Symbol : a symbol (string or function handle) to look up; if omitted, the entire symbol frame
%            will be returned as a struct
%
%   Default : the default value, if the symbol does not exist (note: you may also omit the default
%             if you prefer to get an error in case the symbol does not exist)
%
%   Context : Optionally the execution context (stack) at this point, if known.
%             Can be obtained via try/catch.
%
% Out:
%   Data : the data associated with the symbol, or the symbol itself, if the lookup failed
%
% Examples:
%   % resolve the value of the symbol 'test' against the current dynamic scope
%   value = hlp_resolve('test')
%
% See also:
%   hlp_scope
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-03

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

global tracking;

% get stack entries
if nargin < 3
    stack = dbstack;
    entries = {stack.name};
else
    entries = {context.stack.name};
end

% check for all entries whether they refer to a stack frame...
marker = 'make_func/@(f,a,frame__';
frames = find(strncmp(entries,marker,23));

% go through all stack frames...
symbol = char(x);
for k=frames
    frameid = entries{k}(24:end-14);
    % and check if the symbol is present there, then look up
    if isfield(tracking.stack.frames.(frameid),symbol)
        outargs = tracking.stack.frames.(frameid).(symbol);
        return;
    end
end

try
    if isfield(tracking.stack.base,symbol)
        % symbol is fetched from the base workspace
        outargs = tracking.stack.base.(symbol);
    else
        % symbol is not in the base workspace: remains unsubstituted
        outargs = default;
    end
catch
    % base workspace doesn't exist yet: create it...
    tracking.stack.base = struct();
    outargs = default;
end
