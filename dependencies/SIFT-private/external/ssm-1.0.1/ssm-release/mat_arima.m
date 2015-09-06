function [Z T R P1 Tmmask Rmmask P1mmask] = mat_arima(p, d, q, mean)

%MAT_ARIMA Create base matrices for ARIMA(p, d, q) models.
%   [Z T R P1 Tmmask Rmmask P1mmask] = MAT_ARIMA(p, d, q[, mean])
%       Set mean to true to create a model with mean.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

[Z T R P1 Tmmask Rmmask P1mmask] = mat_arma(p, q, mean);
m       = size(Z, 2);
Z       = [ones(1, d) Z];
T       = [triu(ones(d)) ones(d, 1) zeros(d, m-1); zeros(m, d) T];
R       = [zeros(d, 1); R];
P1      = blkdiag(diag(repmat(Inf, d, 1)), P1);
if ~isempty(Tmmask), Tmmask = [false(d, d+m); false(m, d) Tmmask]; end
if ~isempty(Rmmask), Rmmask = [false(d, 1); Rmmask]; end
P1mmask  = logical(blkdiag(zeros(d), P1mmask));

