function varargout = hlp_getresult(idx,f,varargin)
% Returns the Result-Idx's output of the function, given the supplied arguments.
% Results... = hlp_getresult(Result-Idx, Function, Arguments...)
%
% In:
%   Result-Idx : index of the result (output) of the given function application; can also 
%                be a vector of indices, then re-emitting those outputs as outputs.
%                note: if this is passed as a cell array, the outputs will be wrapped into a cell 
%                      array.
%
%   Function : function to apply to the Arguments
%
%   Arguments... : list of arguments to the function
%
% Out:
%   Results... : one or more outputs, selected according to Result-Idx from the outputs of 
%                Function(Arguments...), and optionally wrapped into a cell array.
%
% Examples:
%   % get the second output of the sort() function, when applied to some data
%   hlp_getresult(2,@sort,[1 4 2 5 2 1 0 6])
%
% See also:
%   hlp_wrapresults
%
%				          		 Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-03-28

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

if isnumeric(idx)
    [tmp{1:max(idx)}] = f(varargin{:});
    varargout = tmp(idx);
elseif iscell(idx)
    [tmp{1:max([idx{:}])}] = f(varargin{:});
    varargout = {tmp([idx{:}])};
else
    error('Unsupported index format (must be either numeric or cell array).');
end