function [m mmask] = mat_dupvar(p, cov, d)

%MAT_DUPVAR Create base matrices for duplicated variance noise.
%   [m mmask] = MAT_DUPVAR(p, cov, d)
%       p is the number of variables.
%       d is the number of duplicates.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

m   = zeros(p*d);
if cov, mmask = logical(kron(ones(p), eye(d)));
else mmask = logical(eye(p*d));
end

