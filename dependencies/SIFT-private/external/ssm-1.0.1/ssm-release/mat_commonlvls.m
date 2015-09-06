function [Z T R a1 P1] = mat_commonlvls(p, A, a)

%MAT_COMMONLVLS Create base matrices for common levels model.
%   [Z T R a1 P1] = MAT_COMMONLVLS(p, A, a)
%       A is a (p-r)*r matrix.
%       a is a p-r column vector.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

r   = size(A, 2);
Z   = [eye(r) zeros(r, p-r); A eye(p-r)];
T   = eye(p);
R   = [eye(r); zeros(p-r, r)];
a1  = [zeros(r, 1); a];
P1  = diag([repmat(Inf, r, 1); zeros(p-r, 1)]);

