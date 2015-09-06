function [fun grad param] = fun_freqspec(p, q, P, nparam, nfreq, freq)

%FUN_FREQSPEC Create update functions for frequency specific models.
%   [fun grad param] = FUN_FREQSPEC(p, q, P, nparam, nfreq, freq)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

s_p             = p;
s_q             = q;
s_P             = P;
s_useb          = (nparam == 4);
s_pmask1        = 1:p;
s_pmask2        = p+1:p+P;
s_pmask3        = p+P+1:p+P+q+nparam-1;
s_pmask4        = p+P+q+nparam;
s_Phimask       = [zeros(1, 11) -1];
s_mask          = true(6, 1);
s_mask(freq)    = false;
s_B             = [ones(6, 1) [-2*cos((1:5)*pi/6)'; 1] zeros(6, 1)];
s_C             = zeros(6, 1);
s_r         = max(s_p+12*s_P, s_q+12+1);
s_T         = [zeros(s_r, 1) [eye(s_r-1); zeros(1, s_r-1)]];
s_R         = [1; zeros(s_r-1, 1)];
s_Tmask     = [[true(s_p+12*s_P, 1); false(s_r-(s_p+12*s_P), 1)] false(s_r, s_r-1)];
s_Rmask     = [false; true(s_q+12, 1); false(s_r-(s_q+12)-1, 1)];
s_I         = eye(s_r^2);
s_T0        = s_I - kron(s_T, s_T);
s_R0        = reshape(s_R*s_R', s_r^2, 1);

    function [Rvec Qvec P1vec] = psi2freqspecma(X)
        X3  = X(s_pmask3);
        % GENSMA
        if s_useb
            if s_q == 0
                Y   = X3./realsqrt(1+X3.^2);
                a   = Y(1);
                c1  = Y(2);
                c2  = Y(3);
            elseif s_q == 1
                Y   = X3./realsqrt(1+X3.^2);
                a   = [(1-Y(2))*Y(1) Y(2)];
                c1  = Y(3);
                c2  = Y(4);
            else
                a   = X3(1:s_q+1);
                Y   = X3(s_q+2:s_q+3)./realsqrt(1+X3(s_q+2:s_q+3).^2);
                c1  = Y(1);
                c2  = Y(2);
            end
            A   = [1 -a];
        else
            if s_q == 0
                Y   = X3./realsqrt(1+X3.^2);
                a   = [];
                c1  = Y(1);
                c2  = Y(2);
            elseif s_q == 1
                Y   = X3./realsqrt(1+X3.^2);
                a   = Y(1);
                c1  = Y(2);
                c2  = Y(3);
            elseif s_q == 2
                Y   = X3./realsqrt(1+X3.^2);
                a   = [(1-Y(2))*Y(1) Y(2)];
                c1  = Y(3);
                c2  = Y(4);
            else
                a   = X3(1:s_q);
                Y   = X3(s_q+1:s_q+2)./realsqrt(1+X3(s_q+1:s_q+2).^2);
                c1  = Y(1);
                c2  = Y(2);
            end
            A   = conv([1 -a], [1 -c1]);
        end
        s_C(s_mask)     = c1;
        s_C(~s_mask)    = c2;
        B               = s_B;
        B(:, 2)         = B(:, 2).*s_C;
        B(1:5, 3)       = s_C(1:5).^2;
        A               = conv(A, B(1, :));
        A               = conv(A, B(2, :));
        A               = conv(A, B(3, :));
        A               = conv(A, B(4, :));
        A               = conv(A, B(5, :));
        A               = conv(A, B(6, 1:2));
        theta           = A(2:s_q+13);
        % Sigma
        zetavar     = exp(2*X(s_pmask4));
        % P1
        s_R(s_Rmask)    = theta;
        P1              = zetavar*linsolve(s_T0, reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Rvec    = theta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

    function [Tvec Rvec Qvec P1vec] = psi2freqspec(X)
        X1  = X(s_pmask1);
        X2  = X(s_pmask2);
        X3  = X(s_pmask3);
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
        % GENSMA
        if s_useb
            if s_q == 0
                Y   = X3./realsqrt(1+X3.^2);
                a   = Y(1);
                c1  = Y(2);
                c2  = Y(3);
            elseif s_q == 1
                Y   = X3./realsqrt(1+X3.^2);
                a   = [(1-Y(2))*Y(1) Y(2)];
                c1  = Y(3);
                c2  = Y(4);
            else
                a   = X3(1:s_q+1);
                Y   = X3(s_q+2:s_q+3)./realsqrt(1+X3(s_q+2:s_q+3).^2);
                c1  = Y(1);
                c2  = Y(2);
            end
            A   = [1 -a];
        else
            if s_q == 0
                Y   = X3./realsqrt(1+X3.^2);
                a   = [];
                c1  = Y(1);
                c2  = Y(2);
            elseif s_q == 1
                Y   = X3./realsqrt(1+X3.^2);
                a   = Y(1);
                c1  = Y(2);
                c2  = Y(3);
            elseif s_q == 2
                Y   = X3./realsqrt(1+X3.^2);
                a   = [(1-Y(2))*Y(1) Y(2)];
                c1  = Y(3);
                c2  = Y(4);
            else
                a   = X3(1:s_q);
                Y   = X3(s_q+1:s_q+2)./realsqrt(1+X3(s_q+1:s_q+2).^2);
                c1  = Y(1);
                c2  = Y(2);
            end
            A   = conv([1 -a], [1 -c1]);
        end
        s_C(s_mask)     = c1;
        s_C(~s_mask)    = c2;
        B               = s_B;
        B(:, 2)         = B(:, 2).*s_C;
        B(1:5, 3)       = s_C(1:5).^2;
        A               = conv(A, B(1, :));
        A               = conv(A, B(2, :));
        A               = conv(A, B(3, :));
        A               = conv(A, B(4, :));
        A               = conv(A, B(5, :));
        A               = conv(A, B(6, 1:2));
        theta           = A(2:s_q+13);
        zetavar     = exp(2*X(s_pmask4));
        % P1
        s_T(s_Tmask)    = TotalPhi;
        s_R(s_Rmask)    = theta;
        P1              = zetavar*linsolve(s_I-kron(s_T, s_T), reshape(s_R*s_R', s_r^2, 1));
        % Set output
        Tvec    = TotalPhi';
        Rvec    = theta';
        Qvec    = zetavar;
        P1vec   = P1;
    end

if p+12*P == 0, fun = {@psi2freqspecma};
else fun = {@psi2freqspec};
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
if nparam == 4
    group(k)    = q + 1;
    if q == 0, transform{k} = 'arma1';
    elseif q == 1, transform{k} = 'ar2';
    else transform{k} = 'ar>=3'; end
    k           = k + 1;
else
    if q > 0
        group(k)    = q;
        if q == 1, transform{k} = 'arma1';
        elseif q == 2, transform{k} = 'ar2';
        else transform{k} = 'ar>=3'; end
        k           = k + 1;
    end
end
group(k:k+2)        = [1 1 1];
[transform{k:k+2}]  = deal('arma1', 'arma1', '1/2 log');
psiname = cell(1, p+P+q+nparam);
for i = 1:p, psiname{i}         = ['phi' int2str(i)]; end
for i = 1:P, psiname{p+i}       = ['Phi' int2str(i)]; end
for i = 1:q+nparam-3, psiname{p+P+i} = ['a' int2str(i)]; end
psiname{p+P+q+nparam-2} = 'c1';
psiname{p+P+q+nparam-1} = 'c2';
psiname{p+P+q+nparam}   = 'zeta var';
param   = ssparam(psiname, transform, group);
end

