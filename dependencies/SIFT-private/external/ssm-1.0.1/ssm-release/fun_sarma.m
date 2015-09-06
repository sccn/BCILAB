function [fun grad param] = fun_sarma(p, q, P, Q, s)

%FUN_SARMA Create update functions for SARMA(p, q, P, Q)s.
%   [fun grad param] = FUN_SARMA(p, q, P, Q, s)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_p         = p;
s_q         = q;
s_P         = P;
s_Q         = Q;
s_pmask1    = 1:p;
s_pmask2    = p+1:p+P;
s_pmask3    = p+P+1:p+P+q;
s_pmask4    = p+P+q+1:p+P+q+Q;
s_pmask5    = p+P+q+Q+1;
s_Phimask   = [zeros(1, s-1) -1];
s_Thetamask = [zeros(1, s-1) 1];
s_r         = max(s_p+s_P*s, s_q+s_Q*s+1);
s_T         = [zeros(s_r, 1) [eye(s_r-1); zeros(1, s_r-1)]];
s_R         = [1; zeros(s_r-1, 1)];
s_Tmask     = [[true(s_p+s_P*s, 1); false(s_r-(s_p+s_P*s), 1)] false(s_r, s_r-1)];
s_Rmask     = [false; true(s_q+s_Q*s, 1); false(s_r-(s_q+s_Q*s)-1, 1)];
s_I         = eye(s_r^2);
s_T0        = s_I - kron(s_T, s_T);
s_R0        = reshape(s_R*s_R', s_r^2, 1);

    function [Qvec P1vec] = psi2sarma00(X)
        zetavar = exp(2*X);
        Qvec    = zetavar;
        P1vec   = zetavar;
    end

    function [Tvec Qvec P1vec] = psi2sar(X)
        X1  = X(s_pmask1);
        X2  = X(s_pmask2);
        if s_p == 1, phi = X1/realsqrt(1 + X1^2);
        elseif s_p == 2
            Y       = X1./realsqrt(1 + X1.^2);
            phi     = [(1-Y(2))*Y(1) Y(2)];
        else phi = X1; end
        if s_P == 1, Phi = X2/realsqrt(1 + X2^2);
        elseif s_P == 2
            Y       = X2./realsqrt(1 + X2.^2);
            Phi     = [(1-Y(2))*Y(1) Y(2)];
        else Phi = X2; end
        TotalPhi    = conv([1 -phi], [1 kron(Phi, s_Phimask)]);
        TotalPhi    = -TotalPhi(2:end);
        zetavar     = exp(2*X(s_pmask5));
        % P1
        s_T(s_Tmask)    = TotalPhi;
        P1              = zetavar*linsolve(s_I-kron(s_T, s_T), s_R0);
        % Set output
        Tvec    = TotalPhi';
        Qvec    = zetavar;
        P1vec   = P1;
    end

    function [Rvec Qvec P1vec] = psi2sma(X)
        X3  = X(s_pmask3);
        X4  = X(s_pmask4);
        if s_q == 1, theta = X3/realsqrt(1 + X3^2);
        elseif s_q == 2
            Y       = X3./realsqrt(1 + X3.^2);
            theta   = [(1+Y(2))*Y(1) Y(2)];
        else theta = X3; end
        if s_Q == 1, Theta = X4/realsqrt(1 + X4^2);
        elseif s_Q == 2
            Y       = X4./realsqrt(1 + X4.^2);
            Theta   = [(1+Y(2))*Y(1) Y(2)];
        else Theta = X4; end
        TotalTheta  = conv([1 theta], [1 kron(Theta, s_Thetamask)]);
        TotalTheta  = TotalTheta(2:end);
        % Sigma
        zetavar     = exp(2*X(s_pmask5));
        % P1
        s_R(s_Rmask)    = TotalTheta;
        P1              = zetavar*linsolve(s_T0, reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Rvec    = TotalTheta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

    function [Tvec Rvec Qvec P1vec] = psi2sarma(X)
        X1  = X(s_pmask1);
        X2  = X(s_pmask2);
        X3  = X(s_pmask3);
        X4  = X(s_pmask4);
        if s_p == 1, phi = X1/realsqrt(1 + X1^2);
        elseif s_p == 2
            Y       = X1./realsqrt(1 + X1.^2);
            phi     = [(1-Y(2))*Y(1) Y(2)];
        else phi = X1; end
        if s_P == 1, Phi = X2/realsqrt(1 + X2^2);
        elseif s_P == 2
            Y       = X2./realsqrt(1 + X2.^2);
            Phi     = [(1-Y(2))*Y(1) Y(2)];
        else Phi = X2; end
        TotalPhi    = conv([1 -phi], [1 kron(Phi, s_Phimask)]);
        TotalPhi    = -TotalPhi(2:end);
        if s_q == 1, theta = X3/realsqrt(1 + X3^2);
        elseif s_q == 2
            Y       = X3./realsqrt(1 + X3.^2);
            theta   = [(1+Y(2))*Y(1) Y(2)];
        else theta = X3; end
        if s_Q == 1, Theta = X4/realsqrt(1 + X4^2);
        elseif s_Q == 2
            Y       = X4./realsqrt(1 + X4.^2);
            Theta   = [(1+Y(2))*Y(1) Y(2)];
        else Theta = X4; end
        TotalTheta  = conv([1 theta], [1 kron(Theta, s_Thetamask)]);
        TotalTheta  = TotalTheta(2:end);
        zetavar     = exp(2*X(s_pmask5));
        % P1
        s_T(s_Tmask)    = TotalPhi;
        s_R(s_Rmask)    = TotalTheta;
        P1              = zetavar*linsolve(s_I-kron(s_T, s_T), reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Tvec    = TotalPhi';
        Rvec    = TotalTheta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

if p+P*s == 0 && q+Q*s == 0, fun = {@psi2sarma00};
elseif p+P*s == 0, fun = {@psi2sma};
elseif q+Q*s == 0, fun = {@psi2sar};
else fun = {@psi2sarma};
end
grad    = {[]};
k       = 1;
if p > 0
    group(k)    = p;
    if p == 1, transform{k} = 'arma1';
    elseif p == 2, transform{k} = 'ar2';
    else transform{k} = 'ar>=3'; end
    k           = k + 1;
end
if P > 0
    group(k)    = P;
    if P == 1, transform{k} = 'arma1';
    elseif P == 2, transform{k} = 'ar2';
    else transform{k} = 'ar>=3'; end
    k           = k + 1;
end
if q > 0
    group(k)    = q;
    if q == 1, transform{k} = 'arma1';
    elseif q == 2, transform{k} = 'ma2';
    else transform{k} = 'ma>=3'; end
    k           = k + 1;
end
if Q > 0
    group(k)    = Q;
    if Q == 1, transform{k} = 'arma1';
    elseif Q == 2, transform{k} = 'ma2';
    else transform{k} = 'ma>=3'; end
    k           = k + 1;
end
group(k)        = 1;
transform{k}    = '1/2 log';
psiname = cell(1, p+P+q+Q+1);
for i = 1:p, psiname{i}         = ['phi' int2str(i)]; end
for i = 1:P, psiname{p+i}       = ['Phi' int2str(i)]; end
for i = 1:q, psiname{p+P+i}     = ['theta' int2str(i)]; end
for i = 1:Q, psiname{p+P+q+i}   = ['Theta' int2str(i)]; end
psiname{p+P+q+Q+1} = 'zeta var';
param   = ssparam(psiname, transform, group);
end

