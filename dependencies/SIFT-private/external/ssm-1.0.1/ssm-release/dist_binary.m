function dist = dist_binary()

%DIST_BINARY Create SSDIST object for binary distribution.
%   dist = DIST_BINARY()

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

dist    = ssdist(0, @matf_binary, @logpf_binary);

function [H y] = matf_binary(y, theta)
H = ((1+exp(theta)).^2)./exp(theta);
y = theta - (1+exp(theta)) + H.*y;

function logp = logpf_binary(y, theta)
logp = y.*theta - reallog(1 + exp(theta));

