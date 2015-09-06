function [T,x_0] = pinit(W,A,b)
%PINIT Utility init.-procedure for "preconditioned" iterative methods.
%
% T = pinit(W,A)
% [T,x_0] = pinit(W,A,b)
%
% Initialization for `preconditioning' of general-form problems.
% Here, W holds a basis for the null space of L.
%
% Determines the matrix T needed in the iterative routines for
% treating regularization problems in general form.
%
% If b is also specified then x_0, the component of the solution in
% the null space of L, is also computed.

% Reference: P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, 07/29/97.

% Initialization.
[n,nu] = size(W);

% Special treatment of square L.
if (nu==0), T = []; x_0 = zeros(n,1); return, end

% Compute T.
S = pinv(A*W);
T = S*A;

% If required, also compute x_0.
if (nargin==3), x_0 = W*(S*b); end