function [Z T R] = mat_mvhs(p, s)

%MAT_MVHS Create base matrices for multivariate Harrison and Stevens seasonal component.
%   [Z T R] = MAT_MVHS(p, s)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

Z       = kron(eye(p), [1 zeros(1, s-1)]);
T       = kron(eye(p), [zeros(s-1, 1) eye(s-1); 1 zeros(1, s-1)]);
R       = eye(p*s);

