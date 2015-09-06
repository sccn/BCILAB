function [m mmask] = mat_var(p, cov)

%MAT_VAR Create base matrices for Gaussian noise.
%   [m mmask] = MAT_VAR(p, cov)
%       p is the number of variables.
%       cov is a logical scalar or matrix with p rows that specifies the
%           covariance structure. If cov is a matrix, then every column vector
%           selects a covariant subset. If cov is a scalar, then cov = true is
%           equivalent to cov = true(p, 1), which means all variables are
%           correlated, and cov = false is equivalent to cov = eye(p), then
%           all variables are independent.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if isscalar(cov)
    if cov, cov = true(p, 1);
    else cov = logical(eye(p));
    end
end
m       = zeros(p);
mmask   = false(p);
for i = cov, mmask(i, i) = true; end

