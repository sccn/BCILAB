function [Z Zdmmask Zdvec T R] = mat_mvreg(p, x, dep)

%MAT_MVREG Create base matrices for multivariate regression component.
%   [Z Zdmmask Zdvec T R] = MAT_MVREG(p, x, dep)
%       p is the number of observation variables.
%       dep is a p*size(x, 1) logical matrix which specifies the dependence of
%           each observation with each regression variable.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if ~isempty(dep)
    m       = nnz(dep);
    mask    = cell(1, p);
    Zdvec   = cell(1, p);
    for i = 1 : p
        mask{i}     = ones(1, nnz(dep(i, :)));
        Zdvec{i}    = x(dep(i, :), :);
    end
    Zdmmask = logical(blkdiag(mask{:}));
    Zdvec   = vertcat(Zdvec{:});
else
    m       = p*size(x, 1);
    Zdmmask = logical(kron(eye(p), ones(1, size(x, 1))));
    Zdvec   = repmat(x, p, 1);
end
Z   = zeros(p, m);
T   = eye(m);
R   = zeros(m, 0);

