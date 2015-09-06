function dist = dist_negbinomial(k)

%DIST_NEGBINOMIAL Create SSDIST object for negative binomial distribution.
%   dist = DIST_NEGBINOMIAL(k)
%       k is the number of trials at each time point, or a scalar if the
%           number of trials is stationary.

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if nargin < 1, error('ssm:dist_negbinomial:NotEnoughInputs', 'Insufficient input arguments.'); end
if ~isnumeric(k) || size(k, 1) ~= 1, error('ssm:dist_negbinomial:InputError', 'k must be a scalar or row vector.'); end

s_logfactk1 = gammaln(k); % = log(factorial(k-1))

    function [H y] = matf_negbinomial(y, theta)
        H   = ((1-exp(theta)).^2)./(k.*exp(theta));
        y   = theta - (exp(-theta)-1) + H.*y;
    end

    function logp = logpf_negbinomial(y, theta)
        logp    = y.*theta - k.*(theta-reallog(1+exp(theta))) + gammaln(y) - gammaln(y-k+1) - s_logfactk1;
    end

dist    = ssdist(0, @matf_negbinomial, @logpf_negbinomial, size(k, 2));
end

