function [A,b,x] = wing(n,t1,t2)
%WING Test problem with a discontinuous solution.
%
% [A,b,x] = wing(n,t1,t2)
%
% Discretization of a first kind Fredholm integral eqaution with
% kernel K and right-hand side g given by
%    K(s,t) = t*exp(-s*t^2)                       0 < s,t < 1
%    g(s)   = (exp(-s*t1^2) - exp(-s*t2^2)/(2*s)  0 < s   < 1
% and with the solution f given by
%    f(t) = | 1  for  t1 < t < t2
%           | 0  elsewhere.
%
% Here, t1 and t2 are constants satisfying t1 < t2.  If they are
% not speficied, the values t1 = 1/3 and t2 = 2/3 are used.

% Reference: G. M. Wing, "A Primer on Integral Equations of the
% First Kind", SIAM, 1991; p. 109.

% Discretized by Galerkin method with orthonormal box functions;
% both integrations are done by the midpoint rule.

% Per Christian Hansen, IMM, 09/17/92.

% Initialization.
if (nargin==1)
  t1 = 1/3; t2 = 2/3;
else
  if (t1 > t2), error('t1 must be smaller than t2'), end
end
A = zeros(n,n); h = 1/n;

% Set up matrix.
sti = ((1:n)-0.5)*h;
for i=1:n
  A(i,:) = h*sti.*exp(-sti(i)*sti.^2);
end

% Set up right-hand side.
if (nargout > 1)
  b = sqrt(h)*0.5*(exp(-sti*t1^2)' - exp(-sti*t2^2)')./sti';
end

% Set up solution.
if (nargout==3)
  I = find(t1 < sti & sti < t2);
  x = zeros(n,1); x(I) = sqrt(h)*ones(length(I),1);
end