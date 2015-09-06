function [fun grad param] = fun_oneoverf(m)

%FUN_ONEOVERF Create update functions for one-over-f noise.
%   [fun grad param] = FUN_ONEOVERF(m)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_m     = m;
s_mseq0 = 0 : m-1;
s_mseq1 = 1 : m;
s_T     = [zeros(m, 1) [eye(m-1); zeros(1, m-1)]];
s_R     = [1; zeros(m-1, 1)];
s_R     = reshape(s_R*s_R', m^2, 1);
s_Tmask = [true(m, 1) false(m, m-1)];

    function [Tvec Qvec P1vec] = psi2oneoverf(X)
        alpha   = 1 + X(1)/realsqrt(1+X(1)^2);
        phi     = -cumprod((s_mseq0-alpha/2)./s_mseq1);
        zetavar = exp(2*X(2));
        s_T(s_Tmask)    = phi;
        P1              = zetavar*linsolve(eye(s_m^2)-kron(s_T, s_T), s_R);
        
        % Set output
        Tvec    = phi';
        Qvec    = zetavar;
        P1vec   = P1;
    end

fun     = {@psi2oneoverf};
grad    = {[]};
param   = ssparam({'alpha' 'zeta variance'}, {'alpha' '1/2 log'}, [1 1]);
end


