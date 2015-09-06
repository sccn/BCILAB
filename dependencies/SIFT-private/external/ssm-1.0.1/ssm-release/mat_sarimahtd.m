function [Z T R Qmat P1 Tmmask Rmmask Qmmask P1mmask] = mat_sarimahtd(p, d, q, P, D, Q, s)

%MAT_SARIMAHTD Create base matrices for SARIMA with Hillmer-Tiao decomposition.
%   [Z T R Qmat P1 Tmmask Rmmask Qmmask P1mmask] = MAT_SARIMAHTD(p, d, q, P, D, Q, s)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

bAR     = (p+P*s > 0);
bMA     = (q-p-d+(Q-P-D)*s > 0);
M       = 2 + bAR + bMA;
Z       = cell(1, M);
T       = cell(1, M);
R       = cell(1, M);
P1      = cell(1, M);
Tmmask  = cell(1, M);
Rmmask  = cell(1, M);
P1mmask = cell(1, M);
[Z{1} T{1} R{1} P1{1} Tmmask{1} Rmmask{1} P1mmask{1}] = mat_arima(0, d+D, d+D, false);
[Z{2} T{2} R{2} P1{2} Tmmask{2} Rmmask{2} P1mmask{2}] = mat_sumarma(0, D*(s-1), D, s, false);
m(1)    = size(Z{1}, 2) + size(Z{2}, 2);
k       = 3;
if bAR
    [Z{k} T{k} R{k} P1{k} Tmmask{k} Rmmask{k} P1mmask{k}] = mat_arima(p+P*s, 0, p+P*s, false);
    k   = k + 1;
end
if bMA
    [Z{k} T{k} R{k} P1{k} Tmmask{k} Rmmask{k} P1mmask{k}] = mat_arima(0, 0, q-p-d+(Q-P-D)*s, false);
    m(2)    = size(Z{k}, 2);
end
Z       = [Z{:}];
T       = blkdiag(T{:});
R       = blkdiag(R{:});
Qmat    = zeros(M);
P1      = blkdiag(P1{:});
if bAR
    Tmmask  = blkdiag(zeros(m(1)), Tmmask{3});
    if bMA, Tmmask = blkdiag(Tmmask, zeros(m(2))); end
else Tmmask  = [];
end
Tmmask   = logical(Tmmask);
Rmmask   = logical(blkdiag(Rmmask{:}));
Qmmask   = logical(eye(M));
P1mmask  = logical(blkdiag(P1mmask{:}));
