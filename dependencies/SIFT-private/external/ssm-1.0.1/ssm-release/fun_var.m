function [fun grad param] = fun_var(p, cov, name)

%FUN_VAR Create update functions for Gaussian noise.
%   [fun grad param] = FUN_VAR(p, cov[, name])
%       p is the dimension of the Gaussian distribution.
%       cov specifies complete covariance if true, or complete independence if
%           false.
%       name is optional parameter name, it can be a cell array of names.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 3, name = '';
elseif ischar(name) && ~isempty(name), name = [name ' '];
else for i = 1 : length(name), name{i} = [name{i} ' ']; end
end

if cov
    %% Complete dependence %%
    s_p         = p;
    s_p2        = p*(p-1)/2;
    s_bmmask    = logical(tril(ones(p)) - eye(p));
    [s_i s_j]   = find(s_bmmask);
    s_C         = zeros(p);
    s_pmask1    = 1 : p;
    s_pmask2    = p+1 : p*(p+1)/2;
    s_G         = zeros(p^2, p*(p+1)/2);
end

    % Complete dependence
    function vec = psi2cov(X)
        %%%% TODO: Find covariance parametrization that always preserve positivity
        X2          = X(s_pmask2);
        Y           = exp(X(s_pmask1));
        Y           = Y'*Y;
        C           = s_C;
        C(s_bmmask) = Y(s_bmmask).*(X2./realsqrt(1+X2.^2))';
        C           = C + C' + diag(diag(Y));
        vec         = C(:);
    end

    function [vec grad] = psi2covgrad(X)
        X2          = X(s_pmask2);
        Y1          = exp(X(s_pmask1));
        Y1          = Y1'*Y1;
        Y2          = Y1(s_bmmask);
        Y3          = realsqrt(1+X2.^2)';
        Y4          = Y2.*(Y3.^-3);
        C           = s_C;
        C(s_bmmask) = Y2.*(X2'./Y3);
        C           = C + C' + diag(diag(Y1));
        vec         = C(:);
        grad        = s_G;
        for a = 1 : s_p
            C2          = s_C;
            C2(a, :)    = C(a, :);
            C2(:, a)    = C(:, a);
            C2(a, a)    = 2*C2(a, a);
            grad(:, a)  = C2(:);
        end
        for a = 1 : s_p2
            C2                      = s_C;
            C2(s_i(a), s_j(a))      = Y4(a);
            C2(s_j(a), s_i(a))      = Y4(a);
            grad(:, s_pmask2(a))    = C2(:);
        end
    end

    % Complete independence
    function vec = psi2var(X)
        vec     = exp(2*X)';
    end

    function [vec grad] = psi2vargrad(X)
        vec     = exp(2*X)';
        grad    = diag(2*vec);
    end

if cov && p > 1
    %% Complete dependence %%
    fun     = {@psi2cov};
    grad    = {@psi2covgrad};
    psiname = cell(1, p*(p+1)/2);
    for i = 1 : p, psiname{i} = [name 'var ' int2str(i)]; end
    k = p + 1;
    for i = 1 : p-1, for j = i+1 : p, psiname{k} = [name 'cov ' int2str(i) ' ' int2str(j)]; k = k+1; end, end
    param   = ssparam(psiname, 'covariance');
else
    %% Complete independence %%
    fun     = {@psi2var};
    grad    = {@psi2vargrad};
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

