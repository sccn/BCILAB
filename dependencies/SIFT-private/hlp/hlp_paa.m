function y = hlp_paa(x,bp)

% Compute a Piecewise Aggregate Approximation (PAA) representation. A PAA
% representation for breakpoints b1,b2, ... ,bn is defined as
% 
% y = [mean(x(b1:b2) , mean(x(b2:b3)) , ..., mean(x(bn-1:bn))]
%
% In other words, the PAA is simply the collection of mean values over all
% intervals defined by pairs of breakpoints.
%
% Inputs:
%
%     x:    [num_vars x num_points] multivariate time series
%           
%     bp:   breakpoints for PAA (column indices of x)
%
% Outputs:
%
%     y:    If x is a matrix, y contains PAA taken across all rows 
%           (breakpoints defined for columns)
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

y = zeros(1,length(bp)-1);
for i=1:length(bp)-1
    y(i)  = mean(mean(x(:,bp(i):bp(i+1)),1));
end

        