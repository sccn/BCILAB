function varargout = mat_arma(p, q, mean)

%MAT_ARMA Create base matrices for ARMA(p, q) models.
%   [Z T R P1 Tmask Rmask P1mask] = MAT_ARMA(p, q[, mean])
%   [T R P1 Tmask Rmask P1mask] = MAT_ARMA(p, q[, mean])
%       Set mean to true to create a model with mean.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

r   = max(p, q+1);
if mean
    Z   = [1 zeros(1, r)];
    T   = [zeros(r, 1) eye(r); zeros(1, r) 1];
    R   = [1; zeros(r, 1)];
    P1  = diag([zeros(1, r) Inf]);
    if p == 0, Tmask = []; else Tmask = [[true(p, 1); false(r-p+1, 1)] false(r+1, r)]; end
    if q == 0, Rmask = []; else Rmask = [false; true(q, 1); false(r-q, 1)]; end
    P1mask  = logical(blkdiag(ones(r), 0));
else
    Z   = [1 zeros(1, r-1)];
    T   = [zeros(r, 1) [eye(r-1); zeros(1, r-1)]];
    R   = [1; zeros(r-1, 1)];
    P1  = zeros(r);
    if p == 0, Tmask = []; else Tmask = [[true(p, 1); false(r-p, 1)] false(r, r-1)]; end
    if q == 0, Rmask = []; else Rmask = [false; true(q, 1); false(r-q-1, 1)]; end
    P1mask  = true(r);
end
if nargout == 7, varargout = {Z T R P1 Tmask Rmask P1mask};
else varargout = {T R P1 Tmask Rmask P1mask}; end

