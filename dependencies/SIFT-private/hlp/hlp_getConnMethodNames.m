function [names] = hlp_getConnMethodNames(Conn)
% return the names of all connectivity measures in a SIFT Connectivity
% object
%
% Input:
%       Conn:       SIFT connectivity object
%
% Output
%
%       names:       cell array of names of connectivity measures in Conn
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

if isempty(Conn) || isequal(Conn,1)
    names = {};
    return;
end

names = fieldnames(Conn(1));
names = setdiff_bc(names,{'winCenterTimes','erWinCenterTimes','freqs','mode','options','resampleTrialIdx','dims'});
    
    