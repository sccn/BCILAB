function [Z T R] = mat_mvllt(p)

%MAT_MVLLT Create base matrices for multivariate local level trend model.
%   [Z T R] = MAT_MVLLT(p)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

Z   = kron(eye(p), [1 0]);
T   = kron(eye(p), [1 1; 0 1]);
R   = eye(2*p);

