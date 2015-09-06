function varargout = mat_sarima(p, d, q, P, D, Q, s, mean)

%MAT_SARIMA Create base matrices for SARIMA(p, d, q, P, D, Q)s models.
%    [Z T R P1 Tmmask Rmmask P1mmask] = MAT_SARIMA(p, d, q, P, D, Q, s[, mean])
%    [Z T R P1 Rmmask P1mmask] = MAT_SARIMA(p, d, q, P, D, Q, s[, mean])
%       Set mean to true to create a model with mean.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

TotalD  = d + D*s;
[T R P1 Tmmask Rmmask P1mmask] = mat_arma(p+P*s, q+Q*s, mean);
m       = size(T, 1);
Z       = [ones(1, d+D+1) zeros(1, D*(s-1)+m-1)];
TD      = [triu(ones(d+D)) ones(d+D, 1) zeros(d+D, D*(s-1)+m-1)];
TS      = [repmat(-1, s-1, 1) eye(s-1)];
for i = 1:D, TD = [TD; zeros(s-1, d+D+(i-1)*(s-1)) TS zeros(s-1, (D-i)*(s-1)+m-1)]; end
T       = [TD; zeros(m, TotalD) T];
R       = [zeros(TotalD, 1); R];
P1      = blkdiag(diag(repmat(Inf, TotalD, 1)), P1);
if ~isempty(Tmmask), Tmmask = [false(TotalD, TotalD+m); false(m, TotalD) Tmmask]; end
if ~isempty(Rmmask), Rmmask = [false(TotalD, 1); Rmmask]; end
P1mmask  = logical(blkdiag(zeros(TotalD), P1mmask));
if nargout == 7, varargout   = {Z T R P1 Tmmask Rmmask P1mmask};
else varargout   = {Z T R P1 Rmmask P1mmask};
end

