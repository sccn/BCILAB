function [fun grad param] = fun_arma(p, q)

%FUN_ARMA Create update functions for ARMA(p, q).
%   [fun grad param] = FUN_ARMA(p, q)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_p         = p;
s_q         = q;
s_pmask1    = 1:p;
s_pmask2    = p+1:p+q;
s_pmask3    = p+q+1;
s_r         = max(s_p, s_q+1);
s_T         = [zeros(s_r, 1) [eye(s_r-1); zeros(1, s_r-1)]];
s_R         = [1; zeros(s_r-1, 1)];
s_Tmask     = [[true(s_p, 1); false(s_r-(s_p), 1)] false(s_r, s_r-1)];
s_Rmask     = [false; true(s_q, 1); false(s_r-(s_q)-1, 1)];
s_I         = eye(s_r^2);
s_T0        = s_I - kron(s_T, s_T);
s_R0        = reshape(s_R*s_R', s_r^2, 1);

    function [Qvec P1vec] = psi2arma00(X)
        zetavar = exp(2*X);
        Qvec    = zetavar;
        P1vec   = zetavar;
    end

    function [Tvec Qvec P1vec] = psi2ar(X)
        X1  = X(s_pmask1);
        if s_p == 1, phi = X1/realsqrt(1 + X1^2);
        elseif s_p == 2
            Y       = X1./realsqrt(1 + X1.^2);
            phi     = [(1-Y(2))*Y(1) Y(2)];
        else phi = X1; end
        zetavar = exp(2*X(s_pmask3));
        % P1
        s_T(s_Tmask)    = phi;
        P1              = zetavar*linsolve(s_I - kron(s_T, s_T), s_R0);
        % Set output
        Tvec    = phi';
        Qvec    = zetavar;
        P1vec   = P1;
    end

    function [Rvec Qvec P1vec] = psi2ma(X)
        X2  = X(s_pmask2);
        if s_q == 1, theta = X2/realsqrt(1 + X2^2);
        elseif s_q == 2
            Y       = X2./realsqrt(1 + X2.^2);
            theta   = [(1+Y(2))*Y(1) Y(2)];
        else theta = X2; end
        % Sigma
        zetavar  = exp(2*X(s_pmask3));
        % P1
        s_R(s_Rmask)    = theta;
        P1              = zetavar*linsolve(s_T0, reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Rvec    = theta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

    function [Tvec Rvec Qvec P1vec] = psi2arma(X)
        X1  = X(s_pmask1);
        X2  = X(s_pmask2);
        if s_p == 1, phi = X1/realsqrt(1 + X1^2);
        elseif s_p == 2
            Y       = X1./realsqrt(1 + X1.^2);
            phi     = [(1-Y(2))*Y(1) Y(2)];
        else phi = X1; end
        if s_q == 1, theta = X2/realsqrt(1 + X2^2);
        elseif s_q == 2
            Y       = X2./realsqrt(1 + X2.^2);
            theta   = [(1+Y(2))*Y(1) Y(2)];
        else theta = X2; end
        zetavar  = exp(2*X(s_pmask3));
        % P1
        s_T(s_Tmask)    = phi;
        s_R(s_Rmask)    = theta;
        P1              = zetavar*linsolve(s_I - kron(s_T, s_T), reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Tvec    = phi';
        Rvec    = theta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

if p == 0 && q == 0, fun = {@psi2arma00};
elseif p == 0, fun = {@psi2ma};
elseif q == 0, fun = {@psi2ar};
else fun = {@psi2arma}; end
grad    = {[]};
k       = 1;
if p > 0
    group(k)    = p;
    if p == 1, transform{k} = 'arma1';
    elseif p == 2, transform{k} = 'ar2';
    else transform{k} = 'ar>=3';
    end
    k           = k + 1;
end
if q > 0
    group(k)    = q;
    if q == 1, transform{k} = 'arma1';
    elseif q == 2, transform{k} = 'ma2';
    else transform{k} = 'ma>=3';
    end
    k           = k + 1;
end
group(k)        = 1;
transform{k}    = '1/2 log';
psiname = cell(1, p+q+1);
for i = 1:p, psiname{i} = ['phi' int2str(i)]; end
for i = 1:q, psiname{p+i} = ['theta' int2str(i)]; end
psiname{p+q+1} = 'zeta var';
param   = ssparam(psiname, transform, group);
end

