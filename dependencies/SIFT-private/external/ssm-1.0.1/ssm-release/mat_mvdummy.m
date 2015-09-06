function [Z T R] = mat_mvdummy(p, s, fixed)

%MAT_MVDUMMY Create base matrices for multivariate dummy seasonal component.
%   [Z T R] = MAT_MVDUMMY(p, s[, fixed])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, fixed = false; end

Z   = kron(eye(p), [1 zeros(1, s-2)]);
T   = kron(eye(p), [-ones(1, s-1); eye(s-2) zeros(s-2, 1)]);
if fixed, R = zeros(p*(s-1), 0);
else R = kron(eye(p), [1; zeros(s-2, 1)]);
end

