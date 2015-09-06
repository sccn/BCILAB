function model = ssm_err()

%SSM_ERR Create SSMODEL object for general error distribution noise model.
%   model = SSM_ERR()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

    function distf = psi2err(X)
        %PSI2ERR Update function for general error.
        %   distf = PSI2ERR(X)
        %       X(1) is the variance, X(2) is l.

        epsvar          = exp(2*X(1));
        l               = 1.5 + 0.5*X(2)/realsqrt(1+X(2)^2);
        wl              = 2*realsqrt(gamma(3*l/4))/(l*gamma(l/4)^(3/2));
        cl              = (gamma(3*l/4)/gamma(l/4))^(l/2);
        wlepsvar        = wl/realsqrt(epsvar);
        clepsvar        = cl/epsvar^(l/2);
        N1              = 1/(l*clepsvar);
        logwlepsvar     = reallog(wlepsvar);

        function H = matf_err(eps)
            H       = N1./(abs(eps).^(l-2));
        end

        function logp = logpf_err(eps)
            logp    = logwlepsvar - clepsvar*(abs(eps).^l);
        end

        distf   = {@matf_err @logpf_err};
    end

model   = ssmodel('general error noise', ssdist(1), zeros(1, 0), [], [], [], 'Hng', {@psi2err}, {[]}, ssparam({'error variance' 'error l'}, {'1/2 log' 'l'}, [1 1]));
end

