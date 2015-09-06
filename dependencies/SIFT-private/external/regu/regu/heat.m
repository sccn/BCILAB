function [A,b,x] = heat(n,kappa)
%HEAT Test problem: inverse heat equation.
%
% [A,b,x] = heat(n,kappa)
%
% A first kind Volterra integral equation with [0,1] as
% integration interval.  The kernel is K(s,t) = k(s-t) with
%    k(t) = t^(-3/2)/(2*kappa*sqrt(pi))*exp(-1/(4*kappa^2*t)) .
% Here, kappa controls the ill-conditioning of the matrix:
%    kappa = 5 gives a well-conditioned problem
%    kappa = 1 gives an ill-conditioned problem.
% The default is kappa = 1.
%
% An exact soltuion is constructed, and then the right-hand side
% b is produced as b = A*x.

% Reference: A. S. Carasso, "Determining surface temperatures from interior
% observations", SIAM J. Appl. Math. 42 (1982), 558-574.  See also L. Elden,
% "The numerical solution of a non-characteristic Cauchy problem for a
% parabolic equation"; in P. Deuflhand and E. Hairer (Eds.), "Numerical
% Treatment of Inverse Problems in Differential and Integral Equations",
% Birkhauser, 1983.

% Discretization by means of simple quadrature (midpoint rule).

% Per Christian Hansen, IMM, Sep, 13, 2001.

% Set default kappa.
if (nargin==1), kappa = 1; end

% Initialization.
h = 1/n; t = h/2:h:1;
c = h/(2*kappa*sqrt(pi));
d = 1/(4*kappa^2);

% Compute the matrix A.
k = c*t.^(-1.5).*exp(-d./t);
r = zeros(1,length(t)); r(1) = k(1); A = toeplitz(k,r);

% Compute the vectors x and b.
if (nargout>1)
  x = zeros(n,1);
  for i=1:n/2
    ti = i*20/n;
    if (ti < 2)
      x(i) = 0.75*ti^2/4;
    elseif (ti < 3)
      x(i) = 0.75 + (ti-2)*(3-ti);
    else
      x(i) = 0.75*exp(-(ti-3)*2);
    end
  end
  x(n/2+1:n) = zeros(1,n/2);
  b = A*x;
end