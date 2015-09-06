
function v = hlp_struct2varargin_sift(g)
% Converts a structure to a cell array of ('name',value) pairs.
%
% INPUT:
%       g - a structure
%
% OUTPUT:
%       v - a cell array of ('fieldname', value) pairs for each field in the
%       original struct
% 
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
names = fieldnames(g);
vals = struct2cell(g);
v = cell(1,2*length(names));
v(1:2:end) = names;
v(2:2:end) = vals;

