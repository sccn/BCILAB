function dist = dist_binomial(k)

%DIST_BINOMIAL Create SSDIST object for binomial distribution.
%   dist = DIST_BINOMIAL(k)
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:dist_binomial:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(k) || size(k, 1) ~= 1, error('ssm:dist_binomial:InputError', 'k must be a scalar or row vector.'); end

s_logfactk = gammaln(k+1); % = log(factorial(k))

    function [H y] = matf_binomial(y, theta)
        H   = ((1+exp(theta)).^2)./(k.*exp(theta));
        y   = theta - (1 + exp(theta)) + H.*y;
    end

    function logp = logpf_binomial(y, theta)
        logp    = y.*theta - k.*reallog(1+exp(theta)) + s_logfactk - gammaln(k-y+1) - gammaln(y+1);
    end

dist    = ssdist(0, @matf_binomial, @logpf_binomial, size(k, 2));
end

