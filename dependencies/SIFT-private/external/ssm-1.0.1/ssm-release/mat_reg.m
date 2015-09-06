function [Z Zdmmask Zdvec T R] = mat_reg(x, dyn)

%MAT_REG Create base matrices for regression component.
%   [Z Zdmmask Zdvec T R] = MAT_REG(x[, dyn])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, dyn = false; end

Z       = zeros(1, size(x, 1));
Zdmmask = true(1, size(x, 1));
Zdvec   = x;
T       = eye(size(x, 1));
if dyn, R = eye(size(x, 1));
else R = zeros(size(x, 1), 0);
end

