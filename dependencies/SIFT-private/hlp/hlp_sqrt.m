function Conn = hlp_sqrt(Conn,methods)
%
% Compute the square root of connectivities.
%
% Inputs:
% 
%   Conn        Connectivity structure as computed by est_mvarConnectivity().
%   methods     Cell array of strings denoting which connectivity measures
%               (fields of Conn) to transform.
% Outputs:
%  
%   Conn        Transformed connectivity structure.
%
% See Also: pop_est_mvarConnectivity()
%
% References:
%
% [1] Mullen T (2010) The Source Information Flow Toolbox (SIFT):
%   Theoretical Handbook and User Manual.
%   Available at: http://www.sccn.ucsd.edu/wiki/Sift/
% 
% 
% Author: Tim Mullen 2011, SCCN/INC, UCSD. 
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

    if ~iscell(methods)
        error('methods must be a cell matrix of strings');
    end

    for cnd=1:length(Conn)
        for i=1:length(methods) 
            if ~isfield(Conn(cnd),methods{i})
                if verb, fprintf('Conn(%d).%s does not exist. Skipping...\n',cnd,methods{i}); end
                continue;
            end
            
            Conn(cnd).(methods{i}) = sqrt(Conn(cnd).(methods{i}));
            
        end
    end