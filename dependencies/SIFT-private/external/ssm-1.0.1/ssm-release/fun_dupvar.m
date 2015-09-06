function [fun grad param] = fun_dupvar(p, cov, d, name)

%FUN_DUPVAR Create update functions for duplicated variance noise.
%   [fun grad param] = FUN_DUPVAR(p, cov, d[, name])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, name = '';
elseif ischar(name), name = [name ' '];
else for i = 1 : length(name), name{i} = [name{i} ' ']; end
end

if cov && p > 1
    s_p         = p;
    s_p2        = p*(p-1)/2;
    s_bmmask    = logical(tril(ones(p)) - eye(p));
    [s_i s_j]   = find(s_bmmask);
    s_C         = zeros(p);
    s_pmask1    = 1 : p;
    s_pmask2    = p+1 : p*(p+1)/2;
    s_D         = ones(1, d);
    s_G         = zeros(d*p^2, p*(p+1)/2);
else
    s_d         = d;
end

    function vec = psi2dupvar1(X)
        vec     = repmat(exp(2*X), s_d, 1);
    end

    function [vec grad] = psi2dupvar1grad(X)
        Y       = exp(2*X)';
        vec     = repmat(Y, s_d, 1);
        grad    = repmat(2*Y, s_d, 1);
    end

    function vec = psi2dupvar(X)
        vec     = kron(exp(2*X)', ones(s_d, 1));
    end

    function [vec grad] = psi2dupvargrad(X)
        Y       = exp(2*X)';
        vec     = kron(Y, ones(s_d, 1));
        grad    = kron(diag(2*Y), ones(s_d, 1));
    end

    function vec = psi2dupcov(X)
        %%%% TODO: Find parametrization that always preserve positivity
        X2          = X(s_pmask2);
        Y           = exp(X(s_pmask1));
        Y           = Y'*Y;
        C           = s_C;
        C(s_bmmask) = Y(s_bmmask).*(X2./realsqrt(1+X2.^2))';
        C           = kron(C + C' + diag(diag(Y)), s_D);
        vec         = C(:);
    end

    function [vec grad] = psi2dupcovgrad(X)
        X2          = X(s_pmask2);
        Y1          = exp(X(s_pmask1));
        Y1          = Y1'*Y1;
        Y2          = Y1(s_bmmask);
        Y3          = realsqrt(1+X2.^2)';
        Y4          = Y2.*(Y3.^-3);
        C           = s_C;
        C(s_bmmask) = Y2.*(X2'./Y3);
        C           = C + C' + diag(diag(Y1));
        grad        = s_G;
        for a = 1 : s_p
            C2          = s_C;
            C2(a, :)    = C(a, :);
            C2(:, a)    = C(:, a);
            C2(a, a)    = 2*C2(a, a);
            C2          = kron(C2, s_D);
            grad(:, a)  = C2(:);
        end
        for a = 1 : s_p2
            C2                      = s_C;
            C2(s_i(a), s_j(a))      = Y4(a);
            C2(s_j(a), s_i(a))      = Y4(a);
            C2                      = kron(C2, s_D);
            grad(:, s_pmask2(a))    = C2(:);
        end
        C           = kron(C, s_D);
        vec         = C(:);
    end

if p == 1
    fun     = {@psi2dupvar1};
    grad    = {@psi2dupvar1grad};
    param   = ssparam({[name 'var']}, '1/2 log');
elseif cov
    fun     = {@psi2dupcov};
    grad    = {@psi2dupcovgrad};
    psiname = cell(1, p*(p+1)/2);
    for i = 1 : p, psiname{i} = [name 'var ' int2str(i)]; end
    k = p + 1;
    for i = 1 : p-1, for j = i+1 : p, psiname{k} = [name 'cov ' int2str(i) ' ' int2str(j)]; k = k+1; end, end
    param   = ssparam(psiname, 'covariance');
else % cov == false
    fun     = {@psi2dupvar};
    grad    = {@psi2dupvargrad};
    psiname = cell(1, p);
    if ischar(name)
        if p == 1, psiname{1} = [name 'var'];
        else for i = 1 : p, psiname{i} = [name 'var ' int2str(i)]; end
        end
    else
        for i = 1 : p, psiname{i} = [name{i} 'var']; end
    end
    param   = ssparam(psiname, '1/2 log');
end
end

