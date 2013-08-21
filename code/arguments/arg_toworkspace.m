function arg_toworkspace(args)
% Copy the arguments in the given Struct into the workspace of the calling function. 
% arg_toworkspace(Struct)
%
% In:
%   Struct : an argument structure, as produced by arg_define
%
% Examples:
%   function myfunction(varargin)
%   opts = arg_define(varargin, ...
%               arg(...), ...
%               arg(...));
%
%   arg_toworkspace(opts);
%
% See also:
%   arg_define
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-09-24

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

% place the variables in the caller's workspace
for fn=fieldnames(args)'
    assignin('caller',fn{1},args.(fn{1})); end
