function [dist,cp] = ncpfun(lambda,s,beta,U,dsvd)

% Auxiliary routine for ncp.  PCH, IMM, Dec. 30, 2007.

if (nargin==4)
   f = (lambda^2)./(s.^2 + lambda^2);
else
   f = lambda./(s + lambda);
end
r = U*(f.*beta); m = length(r);
if isreal(beta), q = floor(m/2); else q = m-1; end
D = abs(fft(r)).^2; D = D(2:q+1); v = (1:q)'/q;
cp = cumsum(D)/sum(D);
dist = norm(cp-v);