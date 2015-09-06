function [m mmask] = mat_wvar(p, s)

%MAT_WVAR Create base matrices for W structure variance noise.
%   [m mmask] = MAT_WVAR(p, s)
%       p is the number of variables.
%       s is the size of the W structure.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

m       = zeros(p*s);
mmask   = logical(kron(eye(p), ones(s)));

