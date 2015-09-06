function [H Q Hmmask Qmmask] = mat_homovar(p, cov, q)

%MAT_HOMOVAR Create base matrices for homogeneous variance noise.
%   [H Q Hmmask Qmmask] = MAT_HOMOVAR(p, cov, q)
%       p is the number of variables.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

H   = zeros(p);
Q   = zeros(p*q);
if cov
    Hmmask  = true(p);
    Qmmask  = logical(kron(ones(p), eye(q)));
else
    Hmmask  = logical(eye(p));
    Qmmask  = logical(eye(p*q));
end
