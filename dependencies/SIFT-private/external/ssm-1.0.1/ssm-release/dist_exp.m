function dist = dist_exp()

%DIST_EXP Create SSDIST object for exponential distribution.
%   dist = DIST_EXP()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

dist    = ssdist(0, @matf_exp, @logpf_exp);

function [H y] = matf_exp(y, theta)
H = theta.^2;
y = 2*theta + H.*y;

function logp = logpf_exp(y, theta)
logp = y.*theta + reallog(theta);

