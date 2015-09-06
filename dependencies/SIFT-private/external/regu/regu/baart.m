function [A,b,x] = baart(n)
%BAART Test problem: Fredholm integral equation of the first kind.
%
% [A,b,x] = baart(n)
%
% Discretization of a first-kind Fredholm integral equation with
% kernel K and right-hand side g given by
%    K(s,t) = exp(s*cos(t)) ,  g(s) = 2*sinh(s)/s ,
% and with integration intervals  s in [0,pi/2] ,  t in [0,pi] .
% The solution is given by
%    f(t) = sin(t) .
%
% The order n must be even.

% Reference: M. L. Baart, "The use of auto-correlation for pseudo-
% rank determination in noisy ill-conditioned linear least-squares
% problems", IMA J. Numer. Anal. 2 (1982), 241-247.

% Discretized by the Galerkin method with orthonormal box functions;
% one integration is exact, the other is done by Simpson's rule.

% Per Christian Hansen, IMM, 09/16/92.

% Check input.
if (rem(n,2)~=0), error('The order n must be even'), end

% Generate the matrix.
hs = pi/(2*n); ht = pi/n; c = 1/(3*sqrt(2));
A = zeros(n,n); ihs = (0:n)'*hs; n1 = n+1; nh = n/2;
f3 = exp(ihs(2:n1)) - exp(ihs(1:n));
for j=1:n
  f1 = f3; co2 = cos((j-.5)*ht); co3 = cos(j*ht);
  f2 = (exp(ihs(2:n1)*co2) - exp(ihs(1:n)*co2))/co2;
  if (j==nh)
    f3 = hs*ones(n,1);
  else
    f3 = (exp(ihs(2:n1)*co3) - exp(ihs(1:n)*co3))/co3;
  end
  A(:,j) = c*(f1 + 4*f2 + f3);
end

% Generate the right-hand side.
if (nargout>1)
  si(1:2*n) = (.5:.5:n)'*hs; si = sinh(si)./si;
  b = zeros(n,1);
  b(1) = 1 + 4*si(1) + si(2);
  b(2:n) = si(2:2:2*n-2) + 4*si(3:2:2*n-1) + si(4:2:2*n);
  b = b*sqrt(hs)/3;
end

% Generate the solution.
if (nargout==3)
  x = -diff(cos((0:n)'*ht))/sqrt(ht);
end