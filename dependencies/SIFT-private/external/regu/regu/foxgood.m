function [A,b,x] = foxgood(n)
%FOXGOOD Test problem: severely ill-posed problem.
%
% [A,b,x] = foxgood(n)
%
% This is a model problem which does not satisfy the
% discrete Picard condition for the small singular values.
% The problem was first used by Fox & Goodwin.

% Reference: C. T. H. Baker, "The Numerical Treatment of
% Integral Equations", Clarendon Press, Oxford, 1977; p. 665.

% Discretized by simple quadrature (midpoint rule).

% Per Christian Hansen, IMM, 03/16/93.

% Initialization.
h = 1/n; t = h*((1:n)' - 0.5);

A = h*sqrt((t.^2)*ones(1,n) + ones(n,1)*(t.^2)');
x = t; b = ((1+t.^2).^1.5 - t.^3)/3;
