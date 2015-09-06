function phi_t = stat_ubsplsm_mkspl(knots,N,order,verb)
% Construct orthonormal univariate smoothing splines
%
% Input                    Information
% -------------------------------------------------------------------------
% knots                    Cector containing index positions of knots along
%                          abscissa (1:N)
%
% N                        Length of the time-series to smooth
%
% order                    B-spline order. Usually this is 4
%
% verb                     Verbosity level. 
%                          0 = no output, 1 = text, 2 = graphical   
%
% Output                   Information
% -------------------------------------------------------------------------
% phi_t                    [K x N] matrix containing K basis functions
%                          centered at each of K (augmented) knot
%                          locations. Generally, K=length(knots)+2 (we add
%                          a knot at each endpoint of the sequence)
%
% See Also: stat_bsplSmoother()
%
%
% Author: Tim Mullen and Wes Thompson, 2010-12, SCCN/INC, UCSD.
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


if nargin<3 || isempty(order)
    order = 4;
end
if nargin<4
    verb = 2;
end

if verb==2
    multiWaitbar('Building Splines','Reset','Color',hlp_getNextUniqueColor);
end

K = length(knots)+2;
abscissa=1:N;

% construct spline basis functions
phi_t=spcol(augknt(knots,order),order,abscissa)';

% orthonormalize basis functions
% TODO: in future, this can be replaced by (qr() or orth())
for k = 1:K
    
    if k>1
        for kk = 1:(k-1)
            % orthogonalize basis functions
            phi_t(k,:) = phi_t(k,:) ...
                         - ((phi_t(k,:)*phi_t(kk,:)')/norm(phi_t(kk,:)).^2) ...
                         * phi_t(kk,:);
        end
    end
    
    % normalize basis functions by euclidean length
    phi_t(k,:) = phi_t(k,:)/norm(phi_t(k,:));
    
    if verb==2
        multiWaitbar('Building Splines',k/K);
    end
end

