function [A,b] = ursell(n)
%URSELL Test problem: integral equation wiht no square integrable solution.
%
% [A,b] = ursell(n)
%
% Discretization of a first kind Fredholm integral equation with
% kernel K and right-hand side g given by
%    K(s,t) = 1/(s+t+1) ,  g(s) = 1 ,
% where both integration itervals are [0,1].
%
% Note: this integral equation has NO square integrable solution.

% Reference: F. Ursell, "Introduction to the theory of linear
% integral equations", Chapter 1 in L. M. Delves & J. Walsh (Eds.),
% "Numerical Solution of Integral Equations", Clarendon Press, 1974.

% Discretized by Galerkin method with orthonormal box functions.

% Per Christian Hansen, IMM, 09/16/92.

% Compute the matrix A.
r = zeros(n,1); c = r;
for k = 1:n
  d1 = 1 + (1+k)/n; d2 = 1 + k/n; d3 = 1 + (k-1)/n;
  c(k) = n*(d1*log(d1) + d3*log(d3) - 2*d2*log(d2));
  e1 = 1 + (n+k)/n; e2 = 1 + (n+k-1)/n; e3 = 1 + (n+k-2)/n;
  r(k) = n*(e1*log(e1) + e3*log(e3) - 2*e2*log(e2));
end
A = hankel(c,r);

% Compute the right-hand side b.
b = ones(n,1)/sqrt(n);