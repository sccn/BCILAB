function result = hlp_iscaller(func, level)
% Test whether some function is calling this function at some level(s) of indirection.
% Result = hlp_iscaller(Function, Levels)
%
% It can be specified what levels of nesting outside the function which runs hlp_iscaller are considered.
%
% In:
%   Function : function handle to a function to be tested
%   Levels   : nesting level(s) that shall be tested (default: all)
%              level 1 is the function that invokes hlp_iscaller
%
% Out: 
%   Result   : whether the current code is called by the Caller
%
% Examples:
%   % in a function, test if the calling function, or the caller of the calling function, is named
%   % 'somefunction'
%   hlp_iscaller('somefunction',1:2)
%
% See also:
%   hlp_getcaller
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-14

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

try
    throw; %#ok<LTARG> % fastest way to get an exception
catch context
    if ~exist('level','var')
        level = 2:length(context.stack); end
    level = level+1;
    level(level < 1 | level > length(context.stack)) = [];
    result = any(strcmp(char(func),{context.stack(level).name}));
end
