function [fun grad param] = fun_spline(delta)

%FUN_SPLINE Create update functions for cubic spline smoothing.
%   [fun grad param] = FUN_SPLINE(delta)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_delta  = [delta.^3/3; delta.^2/2; delta.^2/2; delta];

    function dsubvec = psi2deltaQ(X)
        dsubvec = exp(2*X)*s_delta;
    end

fun     = {@psi2deltaQ};
grad    = {[]};
param   = ssparam({'zeta var'}, '1/2 log');
end

