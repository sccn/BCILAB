function flat = hlp_flattenStruct(s,depth,exclude)
%
% 'flatten' a structure. Recursively reassign all fields of any
% substructure to the parent structure. If duplicate fieldnames exist, they
% will be overwritten
%
% Input:
%       s:          the input structure
%       depth:      optional recursion depth {default: Inf}
%       exclude:    a cell array of fieldnames to exclude from flattening
% Output:
%       flat:   the flattened structure
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
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


    if nargin < 3
        exclude = {};
    elseif ischar(exclude)
        exclude = {exclude};
    end
    
    if nargin < 2
        depth = Inf;
    end
    
    if depth <= 0 || isempty(s)
        % we've reached maximum recursion depth
        flat = s;
        return;
    end
    
    fn = fieldnames(s);
    
    for i=1:length(fn)
        fni = fn{i};
        if isstruct(s.(fni)) && ~ismember_bc({fni},exclude)
            substruct = hlp_flattenStruct(s.(fni),depth-1,exclude);     % recursively flatten the child struct
            fnsub = fieldnames(substruct);
            for j=1:length(fnsub);                      % assign all child fields to parent
                fnsubj = fnsub{j};
                flat.(fnsubj) = substruct.(fnsubj);    
            end
        else
            % assign directly
            flat.(fni) = s.(fni);
        end
    end
