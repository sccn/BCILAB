function [fun grad param] = fun_wvar(p, s, name)

%FUN_WVAR Create update functions for W structure variance noise.
%   [fun grad param] = FUN_WVAR(p, s[, name])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 2, name = ''; else name = [name ' ']; end

s_W     = reshape(eye(s) - repmat(1/s, s), s*s, 1);

    function vec = psi2Wvar1(X)
        vec     = exp(2*X)*s_W;
    end

    function [vec grad] = psi2Wvar1grad(X)
        Y       = exp(2*X);
        vec     = Y*s_W;
        grad    = (2*Y)*s_W;
    end

    function vec = psi2Wvar(X)
        vec     = kron(exp(2*X)', s_W);
    end

    function [vec grad] = psi2Wvargrad(X)
        Y       = exp(2*X)';
        vec     = kron(Y, s_W);
        grad    = kron(diag(2*Y), s_W);
    end

if p == 1
    fun     = {@psi2Wvar1};
    grad    = {@psi2Wvar1grad};
    param   = ssparam({[name 'var']}, '1/2 log');
else
    fun     = {@psi2Wvar};
    grad    = {@psi2Wvargrad};
    psiname = cell(1, p);
    if ischar(name), for i = 1 : p, psiname{i} = [name 'var ' int2str(p)]; end
    else for i = 1 : p, psiname{i} = [name{i} ' var']; end, end
    param   = ssparam(psiname, '1/2 log');
end
end
