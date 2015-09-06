function [fun grad param] = fun_homovar(p, cov, q, name)

%FUN_HOMOVAR Create update functions for homogeneous variance noise.
%   [fun grad param] = FUN_HOMOVAR(p, cov, q[, name])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, name = {'' ''};
elseif ischar(name) && ~isempty(name), name = {[name '1 '] [name '2 ']};
else for i = 1 : length(name), name{i} = [name{i} ' ']; end
end

if cov
    %% Complete dependence %%
    s_p         = p;
    s_p2        = p*(p-1)/2;
    s_bmmask    = logical(tril(ones(p)) - eye(p));
    [s_i s_j]   = find(s_bmmask);
    s_C         = zeros(p);
    s_q         = q;
    s_pmask1    = 1 : p;
    s_pmask2    = p+1 : p*(p+1)/2;
    s_pmask3    = p*(p+1)/2 + 1 : p*(p+1)/2 + q;
    s_Qmmask    = logical(kron(ones(p), eye(q)));
    s_HG        = zeros(p*p, p*(p+1)/2 + q);
    s_QG        = zeros(p*p*q, p*(p+1)/2 + q);
else
    s_pmask1    = 1 : p;
    s_pmask2    = p+1 : p+q;
    s_HG        = zeros(p, p+q);
    s_QG        = zeros(p*q, p+q);
end

    function [Hvec Qvec] = psi2homocov(X)
        X2          = X(s_pmask2);
        Y           = exp(X(s_pmask1));
        Y           = Y'*Y;
        C           = s_C;
        C(s_bmmask) = Y(s_bmmask).*(X2./realsqrt(1+X2.^2))';
        C           = C + C' + diag(diag(Y));
        Hvec        = C(:);
        QC          = kron(C, diag(exp(X(s_pmask3))));
        Qvec        = QC(s_Qmmask);
    end

    function [Hvec Qvec Hgrad Qgrad] = psi2homocovgrad(X)
        X2          = X(s_pmask2);
        X3          = diag(exp(X(s_pmask3)));
        Y1          = exp(X(s_pmask1));
        Y1          = Y1'*Y1;
        Y2          = Y1(s_bmmask);
        Y3          = realsqrt(1+X2.^2)';
        Y4          = Y2.*(Y3.^-3);
        C           = s_C;
        C(s_bmmask) = Y2.*(X2'./Y3);
        C           = C + C' + diag(diag(Y1));
        Hvec        = C(:);
        QC          = kron(C, X3);
        Qvec        = QC(s_Qmmask);
        Hgrad       = s_HG;
        Qgrad       = s_QG;
        for a = 1 : s_p
            C2          = s_C;
            C2(a, :)    = C(a, :);
            C2(:, a)    = C(:, a);
            C2(a, a)    = 2*C2(a, a);
            Hgrad(:, a) = C2(:);
            C2          = kron(C2, X3);
            Qgrad(:, a) = C2(s_Qmmask);
        end
        for a = 1 : s_p2
            C2                      = s_C;
            C2(s_i(a), s_j(a))      = Y4(a);
            C2(s_j(a), s_i(a))      = Y4(a);
            Hgrad(:, s_pmask2(a))   = C2(:);
            C2                      = kron(C2, X3);
            Qgrad(:, s_pmask2(a))   = C2(s_Qmmask);
        end
        for a = 1 : s_q
            C2                      = kron(C, diag(X3(:, a)));
            Qgrad(:, s_pmask3(a))   = C2(s_Qmmask);
        end
    end

    function [Hvec Qvec] = psi2homovar(X)
        V       = exp(2*X(s_pmask1))';
        Hvec    = V;
        Qvec    = kron(V, exp(X(s_pmask2))');
    end

    function [Hvec Qvec Hgrad Qgrad] = psi2homovargrad(X)
        X2              = exp(X(s_pmask2))';
        V               = exp(2*X(s_pmask1))';
        Hvec            = V;
        Qvec            = kron(V, X2);
        VD                  = diag(2*V);
        Hgrad               = s_HG;
        Hgrad(:, s_pmask1)  = VD;
        Qgrad               = s_QG;
        Qgrad(:, s_pmask1)  = kron(VD, X2);
        Qgrad(:, s_pmask2)  = kron(V, diag(X2));
    end

if cov && p > 1
    %% Complete dependence %%
    fun     = {@psi2homocov};
    grad    = {@psi2homocovgrad};
    psiname = cell(1, p*(p+1)/2+q);
    for i = 1 : p, psiname{i} = [name{1} 'var ' int2str(i)]; end
    k = p + 1;
    for i = 1 : p-1, for j = i+1 : p, psiname{k} = [name{1} 'cov ' int2str(i) ' ' int2str(j)]; k = k+1; end, end
    k = p*(p+1)/2;
    for i = 1 : q, psiname{k+i} = [name{2} 'q' int2str(i)]; end
    param   = ssparam(psiname, {'covariance' 'log'}, [p*(p+1)/2 q]);
else
    %% Complete independence %%
    fun     = {@psi2homovar};
    grad    = {@psi2homovargrad};
    psiname = cell(1, p+q);
    if length(name) == 2
        if p == 1, psiname{1} = [name{1} 'var'];
        else for i = 1 : p, psiname{i} = [name{1} 'var ' int2str(i)]; end
        end
        if q == 1, psiname{p+1} = [name{2} 'q'];
        else for i = 1 : q, psiname{p+i} = [name{2} 'q' int2str(i)]; end
        end
    else
        for i = 1 : p, psiname{i}   = [name{i} 'var']; end
        for i = 1 : q, psiname{p+i} = [name{p+i} 'q']; end
    end
    param   = ssparam(psiname, {'1/2 log' 'log'}, [p q]);
end
end

