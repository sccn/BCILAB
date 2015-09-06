function G = gcvfun(lambda,s2,beta,delta0,mn,dsvd)

% Auxiliary routine for gcv.  PCH, IMM, Feb. 24, 2008.

% Note: f = 1 - filter-factors.
if (nargin==5)
   f = (lambda^2)./(s2 + lambda^2);
else
   f = lambda./(s2 + lambda);
end
G = (norm(f.*beta)^2 + delta0)/(mn + sum(f))^2;