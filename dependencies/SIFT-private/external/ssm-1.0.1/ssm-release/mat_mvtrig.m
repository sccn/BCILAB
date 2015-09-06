function [Z T R] = mat_mvtrig(p, s, fixed)

%MAT_MVTRIG Create base matrices for multivariate trigonometric seasonal component.
%   [Z T R] = MAT_MVTRIG(p, s[, fixed])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, fixed = false; end

Z   = kron(eye(p), [repmat([1 0], 1, floor((s-1)/2)) ones(1-mod(s, 2))]);
T   = [];
if mod(s, 2) == 0
    for i = 1 : (s/2-1)
        Lambda  = 2*pi*i/s;
        T       = blkdiag(T, [cos(Lambda) sin(Lambda); -sin(Lambda) cos(Lambda)]);
    end
    T   = blkdiag(T, -1);
else % mod(s, 2) == 1
    for i = 1 : ((s-1)/2)
        Lambda  = 2*pi*i/s;
        T       = blkdiag(T, [cos(Lambda) sin(Lambda); -sin(Lambda) cos(Lambda)]);
    end
end
T   = kron(eye(p), T);
if fixed, R = zeros(p*(s-1), 0);
else R = eye(p*(s-1));
end

