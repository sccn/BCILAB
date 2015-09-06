function [fun grad param] = fun_mvcycle(p)

%FUN_MVCYCLE Create update functions for multivariate cycle component.
%   [fun grad param] = FUN_MVCYCLE(p)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_p = p;

    function vec = psi2mvcycle(X)
        Y       = exp(X);
        sinX    = sin(Y);
        cosX    = cos(Y);
        vec     = repmat([cosX; -sinX; sinX; cosX], s_p, 1);
    end

    function [vec grad] = psi2mvcyclegrad(X)
        Y       = exp(X);
        sinX    = sin(Y);
        cosX    = cos(Y);
        vec     = repmat([cosX; -sinX; sinX; cosX], s_p, 1);
        grad    = repmat([-sinX; -cosX; cosX; -sinX], s_p, 1);
    end

fun     = {@psi2mvcycle};
grad    = {@psi2mvcyclegrad};
param   = ssparam({'frequency'}, 'log');
end
