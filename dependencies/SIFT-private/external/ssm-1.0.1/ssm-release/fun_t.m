function [fun grad param] = fun_t(nu)

%FUN_T Create update functions for t-distribution noise.
%   [fun grad param] = FUN_T([nu])
%       nu is fixed if specified to nonzero values.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, nu = 0; end

if nu > 0
    s_nu    = nu;
    s_nu1   = s_nu + 1;
    s_nu12  = s_nu1/2;
    s_N2    = gammaln(s_nu12) - gammaln(s_nu/2);
end

    function distf = psi2t(X)
        %PSI2T Update function for t-distribution w/ fixed degree of freedom.
        %   distf = PSI2T(X)
        %       X(1) is the variance.

        s_epsvarnu2 = exp(2*X(1))*(s_nu-2);
        s_N1        = s_N2 - reallog(s_epsvarnu2)/2;

        function H = matf_t(eps)
            H       = (s_epsvarnu2 + eps.^2)/s_nu1;
        end

        function logp = logpf_t(eps)
            logp    = s_N1 - s_nu12*reallog(1 + eps.^2/s_epsvarnu2);
        end

        distf   = {@matf_t @logpf_t};
    end

    function distf = psi2t0(X)
        %PSI2T0 Update function for t-distribution.
        %   distf = PSI2T0(X)
        %       X(1) is the variance, X(2) is the degree of freedom.

        s_nu        = 2 + exp(X(2));
        s_epsvarnu2 = exp(2*X(1))*(s_nu-2);
        s_nu1       = s_nu + 1;
        s_nu12      = s_nu1/2;
        s_N1        = gammaln(s_nu12) - gammaln(s_nu/2) - reallog(s_epsvarnu2)/2;

        function H = matf_t0(eps)
            H       = (s_epsvarnu2 + eps.^2)/s_nu1;
        end

        function logp = logpf_t0(eps)
            logp    = s_N1 - s_nu12*reallog(1 + eps.^2/s_epsvarnu2);
        end

        distf   = {@matf_t0 @logpf_t0};
    end

if nu > 0
    fun     = {@psi2t};
    grad    = {[]};
    param   = ssparam({'t variance'}, '1/2 log');
else
    fun     = {@psi2t0};
    grad    = {[]};
    param   = ssparam({'t variance' 't degree of freedom'}, {'1/2 log' 'degree of freedom'}, [1 1]);
end
end
