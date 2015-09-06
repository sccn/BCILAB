function [fun grad param] = fun_sarimahtd(p, d, q, P, D, Q, s, gauss)

%FUN_SARIMAHTD Create update functions for SARIMA with Hillmer-Tiao decomposition.
%   [fun grad param] = FUN_SARIMAHTD(p, d, q, P, D, Q, s[, gauss])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 8, gauss = true; end

s_gauss     = gauss;
s_bMA       = (q-p-d+(Q-P-D)*s > 0);
s_p         = p;
s_d         = d;
s_q         = q;
s_P         = P;
s_D         = D;
s_Q         = Q;
s_s         = s;
s_pmask1    = 1:p;
s_pmask2    = p+1:p+P;
s_pmask3    = p+P+1:p+P+q;
s_pmask4    = p+P+q+1:p+P+q+Q;
s_pmask5    = p+P+q+Q+1;
s_Phimask   = [zeros(1, s-1) -1];
s_Thetamask = [zeros(1, s-1) 1];

s_r         = [d+D+1 D*(s-1)+1];
s_T{1}      = [zeros(s_r(1), 1) [eye(s_r(1)-1); zeros(1, s_r(1)-1)]];
s_T{1}      = eye(s_r(1)^2) - kron(s_T{1}, s_T{1});
s_R{1}      = [1; zeros(s_r(1)-1, 1)];
s_Rmask{1}  = [false; true(s_r(1)-1, 1)];

s_T{2}      = [zeros(s_r(2), 1) [eye(s_r(2)-1); zeros(1, s_r(2)-1)]];
s_T{2}      = eye(s_r(2)^2) - kron(s_T{2}, s_T{2});
s_R{2}      = [1; zeros(s_r(2)-1, 1)];
s_Rmask{2}  = [false; true(s_r(2)-1, 1)];

k           = 3;
if p+P > 0
    s_r(k)      = p + P*s + 1;
    s_T{k}      = [zeros(s_r(k), 1) [eye(s_r(k)-1); zeros(1, s_r(k)-1)]];
    s_Tmask     = [[true(s_r(k)-1, 1); false] false(s_r(k), s_r(k)-1)];
    s_R{k}      = [1; zeros(s_r(k)-1, 1)];
    s_Rmask{k}  = [false; true(s_r(k)-1, 1)];
    s_I         = eye(s_r(k)^2);
    k           = k + 1;
end
if s_bMA
    s_r(k)      = q - p - d + (Q - P - D)*s + 1;
    s_T{k}      = [zeros(s_r(k), 1) [eye(s_r(k)-1); zeros(1, s_r(k)-1)]];
    s_T{k}      = eye(s_r(k)^2) - kron(s_T{k}, s_T{k});
    s_R{k}      = [1; zeros(s_r(k)-1, 1)];
    s_Rmask{k}  = [false; true(s_r(k)-1, 1)];
end

    % SMA
    function varargout = psi2htdngma(X)
        X3      = X(s_pmask3);
        X4      = X(s_pmask4);
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
        zetavar         = exp(2*X(s_pmask5));
        TotalTheta      = conv([1 theta], [1 kron(Theta, s_Thetamask)]);
        TotalTheta      = TotalTheta(2:end);
        [theta ksivar]  = htd(s_d + s_D, s_D, s_s, [], TotalTheta, zetavar);

        % P1
        s_R{1}(s_Rmask{1})  = theta{1};
        s_R{2}(s_Rmask{2})  = theta{2};
        P1                  = [ksivar(1)*linsolve(s_T{1}, reshape(s_R{1}*s_R{1}', s_r(1)^2, 1)); ...
                               ksivar(2)*linsolve(s_T{2}, reshape(s_R{2}*s_R{2}', s_r(2)^2, 1))];
        
        if s_bMA
            s_R{3}(s_Rmask{3})  = theta{3};
            P1                  = [P1; ksivar(3)*linsolve(s_T{3}, reshape(s_R{3}*s_R{3}', s_r(3)^2, 1))];
        end

        if s_gauss
            % Set output
            Hvec        = ksivar(end);
            Rvec        = [theta{:}]';
            Qvec        = ksivar(1:end-1)';
            P1vec       = P1;
            varargout   = {Hvec Rvec Qvec P1vec};
        else
            % t-distribution (H):
            % ksivar(end) is the variance, X(end) is the degree of freedom.
            epsvar      = ksivar(end);
            nu          = 2 + exp(X(end));
            epsvarnu2   = epsvar*(nu-2);
            nu1         = nu+1;
            nu12        = nu1/2;
            lambda      = 1/epsvarnu2;
            N1          = gammaln(nu12) - gammaln(nu/2) + 0.5*reallog(lambda);

            % Set output
            Rvec        = [theta{:}]';
            Qvec        = ksivar(1:end-1)';
            P1vec       = P1;
            Hngf        = {@matf_t @logpf_t};
            varargout   = {Rvec Qvec P1vec Hngf};
        end

        function H = matf_t(eps)
            H       = (epsvarnu2 + eps.^2)/nu1;
        end

        function logp = logpf_t(eps)
            logp    = N1 - nu12*reallog(1 + lambda*eps.^2);
        end
    end

    % SARMA
    function varargout = psi2htdng(X)
        X1      = X(s_pmask1);
        X2      = X(s_pmask2);
        X3      = X(s_pmask3);
        X4      = X(s_pmask4);
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
        zetavar         = exp(2*X(s_pmask5));
        TotalPhi        = conv([1 -phi], [1 kron(Phi, s_Phimask)]);
        TotalTheta      = conv([1 theta], [1 kron(Theta, s_Thetamask)]);
        TotalPhi        = TotalPhi(2:end);
        TotalTheta      = TotalTheta(2:end);
        [theta ksivar]  = htd(s_d + s_D, s_D, s_s, TotalPhi, TotalTheta, zetavar);

        % P1
        s_R{1}(s_Rmask{1})  = theta{1};
        s_R{2}(s_Rmask{2})  = theta{2};
        P1                  = [ksivar(1)*linsolve(s_T{1}, reshape(s_R{1}*s_R{1}', s_r(1)^2, 1)); ...
                               ksivar(2)*linsolve(s_T{2}, reshape(s_R{2}*s_R{2}', s_r(2)^2, 1))];
        s_T{3}(s_Tmask)     = TotalPhi;
        s_R{3}(s_Rmask{3})  = theta{3};
        P1                  = [P1; ksivar(3)*linsolve(s_I-kron(s_T{3}, s_T{3}), reshape(s_R{3}*s_R{3}', s_r(3)^2, 1))];

        if s_bMA
            s_R{4}(s_Rmask{4})  = theta{4};
            P1                  = [P1; ksivar(4)*linsolve(s_T{4}, reshape(s_R{4}*s_R{4}', s_r(4)^2, 1))];
        end

        if s_gauss
            % Set output
            Hvec        = ksivar(end);
            Tvec        = TotalPhi';
            Rvec        = [theta{:}]';
            Qvec        = ksivar(1:end-1)';
            P1vec       = P1;
            varargout   = {Hvec Tvec Rvec Qvec P1vec};
        else
            % t-distribution (H):
            % ksivar(end) is the variance, X(end) is the degree of freedom.
            epsvar      = ksivar(end);
            nu          = 2 + exp(X(end));
            epsvarnu2   = epsvar*(nu-2);
            nu1         = nu+1;
            nu12        = nu1/2;
            lambda      = 1/epsvarnu2;
            N1          = gammaln(nu12) - gammaln(nu/2) + 0.5*reallog(lambda);

            % Set output
            Tvec        = TotalPhi';
            Rvec        = [theta{:}]';
            Qvec        = ksivar(1:end-1)';
            P1vec       = P1;
            Hngf        = {@matf_t @logpf_t};
            varargout   = {Tvec Rvec Qvec P1vec Hngf};
        end

        function H = matf_t(eps)
            H       = (epsvarnu2 + eps.^2)/nu1;
        end

        function logp = logpf_t(eps)
            logp    = N1 - nu12*reallog(1 + lambda*eps.^2);
        end
    end

if p+P*s == 0, fun = {@psi2htdngma};
else fun = {@psi2htdng}; end
grad    = {[]};
k   = 1;
if p == 1,      group(k) = 1;   transform{k} = 'arma1'; k = k + 1;
elseif p == 2,  group(k) = 2;   transform{k} = 'ar2';   k = k + 1;
elseif p > 0,   group(k) = p;   transform{k} = 'ar>=3'; k = k + 1;
end
if P == 1,      group(k) = 1;   transform{k} = 'arma1'; k = k + 1;
elseif P == 2,  group(k) = 2;   transform{k} = 'ar2';   k = k + 1;
elseif P > 0,   group(k) = P;   transform{k} = 'ar>=3'; k = k + 1;
end
if q == 1,      group(k) = 1;   transform{k} = 'arma1'; k = k + 1;
elseif q == 2,  group(k) = 2;   transform{k} = 'ma2';   k = k + 1;
elseif q > 0,   group(k) = q;   transform{k} = 'ma>=3'; k = k + 1;
end
if Q == 1,      group(k) = 1;   transform{k} = 'arma1'; k = k + 1;
elseif Q == 2,  group(k) = 2;   transform{k} = 'ma2';   k = k + 1;
elseif Q > 0,   group(k) = Q;   transform{k} = 'ma>=3'; k = k + 1;
end
group(k) = 1;   transform{k} = '1/2 log';

psiname     = cell(1, p+P+q+Q+1);
for i = 1:p, psiname{i}         = ['phi' int2str(i)]; end
for i = 1:P, psiname{p+i}       = ['Phi' int2str(i)]; end
for i = 1:q, psiname{p+P+i}     = ['theta' int2str(i)]; end
for i = 1:Q, psiname{p+P+q+i}   = ['Theta' int2str(i)]; end
psiname{p+P+q+Q+1}  = 'zeta var';

if ~gauss, group(k+1) = 1; transform{k+1} = 'degree of freedom'; psiname{p+P+q+Q+2} = 't degree of freedom'; end

param   = ssparam(psiname, transform, group);
end

