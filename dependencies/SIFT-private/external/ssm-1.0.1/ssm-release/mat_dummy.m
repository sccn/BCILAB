function [Z T R] = mat_dummy(s, fixed)

%MAT_DUMMY Create base matrices for dummy seasonal component.
%   [Z T R] = MAT_DUMMY(s[, fixed])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, fixed = false; end

Z   = [1 zeros(1, s-2)];
T   = [-ones(1, s-1); eye(s-2) zeros(s-2, 1)];
if fixed, R = zeros(s-1, 0);
else R = [1; zeros(s-2, 1)];
end

