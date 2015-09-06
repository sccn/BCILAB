function dist = dist_poisson()

%DIST_POISSON Create SSDIST object for Poisson distribution.
%   dist = DIST_POISSON()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

dist    = ssdist(0, @matf_poisson, @logpf_poisson);

function [H y] = matf_poisson(y, theta)
H = 1./exp(theta);
y = theta - 1 + H.*y;

function logp = logpf_poisson(y, theta)
logp = y.*theta - exp(theta) - gammaln(y+1); % last term = log(factorial(y))

