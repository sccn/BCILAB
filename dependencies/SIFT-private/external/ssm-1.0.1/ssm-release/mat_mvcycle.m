function [Z T Tmmask R] = mat_mvcycle(p)

%MAT_MVCYCLE Create base matrices for multivariate cycle component.
%   [Z T Tmmask R] = MAT_MVCYCLE(p)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

Z       = kron(eye(p), [1 0]);
T       = zeros(2*p);
Tmmask  = logical(kron(eye(p), [1 1; 1 1]));
R       = eye(2*p);

