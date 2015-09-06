function model = ssm_zmsv()

%SSM_ZMSV Create SSMODEL object for zero-mean stochastic volatility error model.
%   model = SSM_ZMSV()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

    function distf = psi2zmsv(X)
        %PSI2ZMSV Update function for zero-mean stochastic volatility error.
        %   distf = PSI2ZMSV(X)

        sigma22     = 2*exp(X)^2;               % 2*sigma^2
        N1          = reallog(pi*sigma22)/2;    % log(2*pi*sigma^2)

        function [H y] = matf_zmsv(y, theta)
            H   = sigma22*exp(theta)./y.^2;
            y   = theta - H/2 + 1;
        end

        function logp = logpf_zmsv(y, theta)
            logp    = -N1-theta/2-(y.^2).*exp(-theta)/sigma22;
        end

        distf   = {@matf_zmsv @logpf_zmsv};
    end

model   = ssmodel('zero-mean stochastic volatility error', ssdist(0), zeros(1, 0), [], [], [], 'Hng', {@psi2zmsv}, {[]}, ssparam({'sigma'}, 'log'));
end
