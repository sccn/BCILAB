function [Z T R P1 Tmmask Rmmask P1mmask] = mat_sumarma(p, q, D, s, mean)

%MAT_SUMARMA Create base matrices for sum integrated ARMA models.
%   [Z T R P1 Tmmask Rmmask P1mmask] = MAT_SUMARMA(p, q, D, s[, mean])
%       Set mean to true to create a model with mean.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if D == 0
    [Z T R P1 Tmmask Rmmask P1mmask] = mat_arma(p, q, mean);
    return;
end
TotalD  = D*(s-1);
[T R P1 Tmmask Rmmask P1mmask] = mat_arma(p, q, mean);
m       = size(T, 1);
Z       = [-1 1 zeros(1, TotalD+m-2)];
TS      = [repmat(-1, s-1, 1) eye(s-1)];
TS2     = [TS zeros(s-1, TotalD-s+m)];
for i = 2:D, TS2 = [TS2; [zeros(s-1, (i-1)*(s-1)) TS zeros(s-1, (D-i)*(s-1)+m-1)]]; end
T       = [TS2; zeros(m, TotalD) T];
R       = [zeros(TotalD, 1); R];
P1      = blkdiag(diag(repmat(Inf, TotalD, 1)), P1);
if ~isempty(Tmmask), Tmmask = [false(TotalD, TotalD+m); false(m, TotalD) Tmmask]; end
if ~isempty(Rmmask), Rmmask = [false(TotalD, 1); Rmmask]; end
P1mmask  = logical(blkdiag(zeros(TotalD), P1mmask));

