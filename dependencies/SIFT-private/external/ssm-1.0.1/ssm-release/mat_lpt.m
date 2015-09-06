function [Z T R] = mat_lpt(d, stochastic)

%MAT_LPT Create base matrices for local polynomial trend model.
%   [Z T R] = MAT_LPT(d[, stochastic])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, stochastic = true; end

Z   = [1 zeros(1, d)];
T   = triu(ones(d+1));
if stochastic, R = eye(d+1);
else R = [zeros(d, 1); 1];
end

