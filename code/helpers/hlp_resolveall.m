function outargs = hlp_resolveall(context)
% Look up all symbols in the current hlp_scope scope.
% Scope = hlp_resolveall(Context)
%
% In:
%   Context : Optionally the execution context (struct with field .stack) at this point, if known.
%             Can be obtained via dbstack or try/catch.
%
% Out:
%   Scope : a struct with values assigned to symbol names
%
% See also:
%   hlp_resolve
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

if nargin < 1
    stack = dbstack;
    entries = {stack.name};
else
    entries = {context.stack.name};
end

% check for all entries whether they refer to a stack frame...
marker = 'make_func/@(f,a,frame__';
frames = find(strncmp(entries,marker,23));

% start with the base scope
try
    scope = tracking.stack.base;
catch
    scope = struct();
    tracking.stack.base = scope;
end
% and walk up the stack, overriding existing symbols
for k=frames(end:-1:1)
    frame = tracking.stack.frames.(entries{k}(24:end-14));
    for fn=fieldnames(frame)'
        scope.(fn{1}) = frame.(fn{1}); end
end
outargs = scope;
