function model = ssm_mix()

%SSM_MIX Create SSMODEL object for Gaussian mixture distribution noise model.
%   model = SSM_MIX()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

    function distf = psi2mix(X)
        %PSI2MIX Update function for Gaussian mixture distribution.
        %   distf = PSI2MIX(X)
        %       X(1) is the variance, X(2) is the ratio, X(3) is the scale.

        s_epsvar        = exp(2*X(1));
        s_lambda        = 1/realsqrt(1+X(2)^2);
        s_chi           = 1 + exp(X(3));
        s_epsvar2       = 2*s_epsvar;
        s_chiepsvar     = s_chi*s_epsvar;
        s_N1            = s_lambda/realsqrt(pi*s_epsvar2);
        s_N2            = (1-s_lambda)/realsqrt(2*pi*s_chiepsvar);
        s_N1epsvar      = s_N1/s_epsvar;
        s_N2chiepsvar   = s_N2/s_chiepsvar;

        function H = matf_mix(eps)
            eps2        = eps.^2/s_epsvar2;
            expeps2     = exp(-eps2);
            expeps2chi  = exp(-eps2/s_chi);
            H           = (s_N1*expeps2 + s_N2*expeps2chi)/(s_N1epsvar*expeps2 + s_N2chiepsvar*expeps2chi);
        end

        function logp = logpf_mix(eps)
            eps2    = eps.^2/s_epsvar2;
            logp    = reallog(s_N1*exp(-eps2) + s_N2*exp(-eps2/s_chi));
        end

        distf   = {@matf_mix @logpf_mix};
    end

model   = ssmodel('Gaussian mixture noise', ssdist(1), zeros(1, 0), [], [], [], 'Hng', {@psi2mix}, {[]}, ssparam({'mixture variance' 'mixture ratio' 'mixture scale'}, {'1/2 log' 'ratio' 'scale'}, [1 1 1]));
end

