function [ssub sin] = hlp_splitstruct(sin,subset)
% Split a structure into two new structs, each containing a subset of the
% original fields. The union of the two structs recovers the original
% structure.
%
% [ssub srem] = hlp_splitstruct(sin,subset);
%
% In: 
%   sin:    A structure
%
%   subset: A cell array of fieldnames to extract from sin
%
% Out:
%   ssub:   A new structure containing fields from subset
%   sin:    Input structure, with subset removed
%
% See Also: 
% 
% 
% Author: Tim Mullen 2009, SCCN/INC, UCSD. 
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

   for i = 1:length(subset)
       if isfield(sin,subset{i})
           ssub.(subset{i}) = sin.(subset{i});
           sin = rmfield(sin,subset{i});
       end
   end