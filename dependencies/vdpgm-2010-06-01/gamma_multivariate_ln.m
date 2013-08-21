function val = gamma_multivariate_ln(x,p);
% function val = gamma_multivariate_ln(x,p);
%
% x: array(1,K)
% p: scalor
%
% x must be more than (p-1)/2
% x should be more than p/2
%
% Gamma_p(x) = pi^(p(p-1)/4) prod_(j=1)^p Gamma(x+(1-j)/2)
% log Gamma_p(x) = p(p-1)/4 log pi + sum_(j=1)^p log Gamma(x+(1-j)/2)

K = length(x);
gammaln_val = gammaln(repmat(x,p,1)+0.5*(1-repmat([1:p]',1,K))); % p by K
val = p*(p-1)*0.25 * log(pi) + sum(gammaln_val,1);


% Local Variables: ***
% mode: matlab ***
% End: ***
