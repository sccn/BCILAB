function [Z T R] = mat_hs(s)

%MAT_HS Create base matrices for Harrison and Stevens seasonal component.
%   [Z T R] = MAT_HS(s)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

Z   = [1 zeros(1, s-1)];
T   = [zeros(s-1, 1) eye(s-1); 1 zeros(1, s-1)];
R   = eye(s);

