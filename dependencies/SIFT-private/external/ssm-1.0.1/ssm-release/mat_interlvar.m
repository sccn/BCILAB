function [m mmask] = mat_interlvar(p, q, cov)

%MAT_INTERLVAR Create base matrices for q-interleaved variance noise.
%   [m mmask] = MAT_INTERLVAR(p, q, cov)
%       p is the number of variables.
%       q is the number of variances affecting each variable.
%       cov is a logical vector that specifies whether each q variances covary
%           across variables.
%       The variances affecting any given variable is always assumed to be
%           independent.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if isscalar(cov), cov = repmat(cov, 1, q); end

bmmask  = logical(tril(ones(p)));
m       = zeros(p*q);
mmask   = false(p*q);
for i = 1 : q
    emask       = false(q);
    emask(i, i) = true;
    if cov(i), mmask = mmask | logical(kron(bmmask, emask));
    else mmask = mmask | logical(kron(eye(p), emask));
    end
end
mmask   = mmask | mmask';

