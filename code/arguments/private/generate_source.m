function source = generate_source(fmt,source)
% Generate a Source function (a function that uses arg_define to parse arguments)
% Source = generate_source(Format,Source)
%
% This is a helper function to handle the Source and Format arguments of the expand_argsub*
% functions.
%
% In:
%   Format : a format specifier, with the same meaning as in arg_define
%
%   Source : either a cell array of argument declarations (as in arg_define) or a function
%            handle that uses arg_define to parse its arguments
%   
% Out:
%   Source : a function that uses arg_define to parse its arguments
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-19

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

if iscell(source)
    % source is a cell array instead of a function: wrap it by an anonymous arg_define-using function
    source = @(varargin) arg_define(fmt,varargin,source{:});
elseif isa(fmt,'function_handle')
    % source is already a function and fmt is a function, too (a pre-parser)
    source = @(varargin) source(fmt(varargin));
elseif ~isequal(fmt,[0 Inf])
    % source is a function and format attempts to redefine its number of positional arguments: we do not handle this case
    error('When Source is a function then Format must be a function or [0 Inf].');
end
