function [A,b] = parallax(n)
%PARALLAX Stellar parallax problem with 28 fixed, real observations.
%
% [A,b] = parallax(n)
%
% Stellar parallax problem with 28 fixed, real observations.
%
% The underlying problem is a Fredholm integral equation of the
% first kind with kernel
%    K(s,t) = (1/sigma*sqrt(2*pi))*exp(-0.5*((s-t)/sigma)^2) ,
% and it is discretized by means of a Galerkin method with n
% orthonormal basis functions.  The right-hand side consists of
% a measured distribution function of stellar parallaxes, and its
% length is fixed, m = 26.  The exact solution, which represents
% the true distribution of stellar parallaxes, in not known.

% Reference: W. M. Smart, "Stellar Dynamics", Cambridge University Press,
% 1938; p. 30.

% Discretized by Galerkin method with orthonormal box functions;
% 2-D integration is done by means of the computational molecule:
%       1   4   1
%       4  16   4
%       1   4   1

% Per Christian Hansen, IMM, 09/16/92.

% Initialization.
a = 0; b = 0.1; m = 26; sigma = 0.014234;
hs = 0.130/m; hx = (b-a)/n; hsh = hs/2; hxh = hx/2;
ss = (-0.03 + [0:m-1]'*hs)*ones(1,n);
xx = ones(m,1)*(a + [0:n-1]*hx);

% Set up the matrix.
A =     16*exp(-0.5*((ss+hsh - xx-hxh)/sigma).^2);
A = A + 4*(exp(-0.5*((ss+hsh - xx    )/sigma).^2) + ...
           exp(-0.5*((ss+hsh - xx-hx )/sigma).^2) + ...
           exp(-0.5*((ss     - xx-hxh)/sigma).^2) + ...
           exp(-0.5*((ss+hs  - xx-hxh)/sigma).^2));
A = A +   (exp(-0.5*((ss     - xx    )/sigma).^2) + ...
           exp(-0.5*((ss+hs  - xx    )/sigma).^2) + ...
           exp(-0.5*((ss     - xx-hx )/sigma).^2) + ...
           exp(-0.5*((ss+hs  - xx-hx )/sigma).^2));
A = sqrt(hs*hx)/(36*sigma*sqrt(2*pi))*A;

% Set up the normalized right-hand side.
b = [3;7;7;17;27;39;46;51;56;50;43;45;43;32;33;29;...
     21;12;17;13;15;12;6;6;5;5]/(sqrt(hs)*640);