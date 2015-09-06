function var = hlp_mergeVarargin(varargin)
% var = hlp_mergeVarargin(varargin)
% Checks if first argument in varargin is a struct and, if so, converts to
% arglist format (cell array of name-value pairs) and appends any remaining
% varargin arguments. Use in combination with finputcheck() to convert
% mixed-type varargin to structure format
%
% This function will be deprecated in SIFT 1.0-beta
%
% Example:
%   var = hlp_mergeVarargin(varargin);
%   g = finputcheck(var, hlp_getDefaultArglist(myPrefix), myFuncName,'ignore');
% 
% Author: Tim Mullen 2010, SCCN/INC, UCSD. 
% Email:  tim@sccn.ucsd.edu

% This function is part of the Source Information Flow Toolbox (SIFT)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if isempty(varargin)
    var = {};
    return;
end

if isstruct(varargin{1})
    var = hlp_struct2varargin(varargin{1});
    varargin = varargin(2:end);
else
    var = {};
end

if length(varargin)>1
    var = [var,varargin];
end