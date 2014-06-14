function varargout = arg_toworkspace(s)
% Copy the arguments in the given Struct into the workspace of the calling function. 
% Outputs... = arg_toworkspace(Struct)
%
% This function can be used to populate the function's workspace with fields of a struct, and
% optionally allows to assign the first k fields of the struct to explicitly named output variables.
% For all other fields variables with the same name as the field are created implicitly in the
% function's workspace by arg_toworkspace.
%
% Using explicit names and return values allows to work around situations where MATLAB "doesn't
% realize" that a given variable has been written to the workspace and accidentally looks up the
% identifier from the MATLAB path (e.g., a function with the same name). The latter issue is
% explained in online posts about "poofing" variables into the workspace.
%
% In:
%   Struct : an argument structure, as produced by arg_define
%
% Out:
%   Outputs... : optional output arguments that shall receive the first k fields in the Struct
%                all remaining fields are written into the workspace by their name.
%
% Examples:
%   function myfunction(varargin)
%   opts = arg_define(varargin, ...
%               arg('myarg1',...), ...
%               arg('myarg2',...), ...
%               arg('myarg3',...));
%
%   % copy all args into the workspace
%   arg_toworkspace(opts);
%
%   % alternative: accept the first two arguments (myarg1,myarg2) into variables named myfirstarg, 
%   % and mysecondarg, copy the rest into the workspace by name (myarg3)
%   [myfirstarg,mysecondarg] = arg_toworkspace(opts);
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


% assign the first k arguments to the given output variables
if nargout > 0
    fields = struct2cell(s);
    varargout(1:nargout) = fields(1:nargout);
end

% place all remaining variables in the workspace by their name
fnames = fieldnames(s);
for fn=fnames(nargout+1:end)'
    assignin('caller',fn{1},s.(fn{1})); end
