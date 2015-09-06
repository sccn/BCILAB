function [fun grad param] = fun_interlvar(p, q, cov, name)

%FUN_INTERLVAR Create update functions for q-interleaved variance noise.
%   [fun grad param] = FUN_INTERLVAR(p, q, cov[, name])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 4, name = '';
elseif ischar(name)
    temp    = cell(1, q);
    for i = 1 : q, temp{i} = [name int2str(i) ' ']; end
    name    = temp;
else for i = 1 : q, name{i} = [name{i} ' ']; end
end
if isscalar(cov), cov = repmat(cov, 1, q); end

    s_p         = p;
    s_p2        = p*(p-1)/2;
    s_bmmask    = logical(tril(ones(p)) - eye(p));
    [s_i s_j]   = find(s_bmmask);
    s_C         = zeros(p);
s_q             = q;
s_cov           = cov;
s_mmask         = false(p*q);
s_emmask        = cell(1, q);
[s_emmask{:}]   = deal(false(q));
s_pmask1        = cell(1, q);
s_pmask2        = cell(1, q);
k               = 0;
for i = 1 : q
    s_emmask{i}(i, i)   = true;
    s_pmask1{i}         = k+1 : k+p;
    if cov(i)
        s_pmask2{i} = k+p+1 : k+p*(p+1)/2;
        k           = k + p*(p+1)/2;
        s_emmask{i} = logical(kron(s_bmmask | eye(p), s_emmask{i}));
        s_emmask{i} = s_emmask{i} | s_emmask{i}';
    else
        s_pmask2{i} = [];
        k           = k + p;
        s_emmask{i} = logical(kron(eye(p), s_emmask{i}));
    end
    s_mmask     = s_mmask | s_emmask{i};
end
s_vmask         = cell(1, q);
for i = 1 : q
    s_vmask{i}  = s_emmask{i}(s_mmask);
end
s_V             = zeros(p*q + nnz(cov)*p*(p-1), 1);
s_G             = zeros(p*q + nnz(cov)*p*(p-1), p*q + nnz(cov)*p*(p-1)/2);

    function vec = psi2interlvar(X)
        vec     = s_V;
        for i = 1 : s_q
            X1  = X(s_pmask1{i});
            if s_cov(i)
                X2              = X(s_pmask2{i})';
                Y               = exp(X1);
                Y               = Y'*Y;
                C               = s_C;
                C(s_bmmask)     = Y(s_bmmask).*(X2./realsqrt(1+X2.^2));
                C               = C + C' + diag(diag(Y));
                vec(s_vmask{i}) = C(:);
            else
                vec(s_vmask{i}) = exp(2*X1);
            end
        end
    end

    function [vec grad] = psi2interlvargrad(X)
        vec     = s_V;
        grad    = s_G;
        for i = 1 : s_q
            X1  = X(s_pmask1{i});
            if s_cov(i)
                X2              = X(s_pmask2{i})';
                Y1              = exp(X1);
                Y1              = Y1'*Y1;
                Y2              = Y1(s_bmmask);
                Y3              = realsqrt(1+X2.^2);
                Y4              = Y2.*(Y3.^-3);
                C               = s_C;
                C(s_bmmask)     = Y2.*(X2./Y3);
                C               = C + C' + diag(diag(Y1));
                vec(s_vmask{i}) = C(:);
                for a = 1 : s_p
                    C2                                  = s_C;
                    C2(a, :)                            = C(a, :);
                    C2(:, a)                            = C(:, a);
                    C2(a, a)                            = 2*C2(a, a);
                    grad(s_vmask{i}, s_pmask1{i}(a))    = C2(:);
                end
                for a = 1 : s_p2
                    C2                                  = s_C;
                    C2(s_i(a), s_j(a))                  = Y4(a);
                    C2(s_j(a), s_i(a))                  = Y4(a);
                    grad(s_vmask{i}, s_pmask2{i}(a))    = C2(:);
                end
            else
                V                               = exp(2*X1);
                vec(s_vmask{i})                 = V;
                grad(s_vmask{i}, s_pmask1{i})   = diag(2*V);
            end
        end
    end

if p == 1
    [fun grad param] = fun_var(q, false, name);
else
    fun         = {@psi2interlvar};
    grad        = {@psi2interlvargrad};
    psiname     = {};
    group       = zeros(1, q);
    transform   = cell(1, q);
    for i = 1 : q
        if cov(i)
            psiname1    = cell(1, p*(p+1)/2);
            for j = 1 : p, psiname1{j} = [name{i} 'var ' int2str(j)]; end
            l   = p + 1;
            for j = 1 : p-1
                for k = j+1 : p
                    psiname1{l} = [name{i} 'cov ' int2str(j) ' ' int2str(k)];
                    l = l + 1;
                end
            end
            group(i)        = p*(p+1)/2;
            transform{i}    = 'covariance';
        else
            psiname1        = cell(1, p);
            for j = 1 : p, psiname1{j} = [name{i} 'var ' int2str(j)]; end
            group(i)        = p;
            transform{i}    = '1/2 log';
        end
        psiname = [psiname psiname1];
    end
    param       = ssparam(psiname, transform, group);
end
end

