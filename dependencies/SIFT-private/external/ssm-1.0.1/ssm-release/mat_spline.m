function [Z T Tdmask Tdvec R] = mat_spline(delta)

%MAT_SPLINE Create base matrices for cubic spline smoothing.
%   [Z T Tdmask Tdvec R] = MAT_SPLINE(delta)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

Z       = [1 0];
T       = [1 0; 0 1];
Tdmask  = logical([0 1; 0 0]);
Tdvec   = delta;
R       = [1 0; 0 1];
