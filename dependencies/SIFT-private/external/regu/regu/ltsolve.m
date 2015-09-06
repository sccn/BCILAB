function x = ltsolve(L,y,W,T)
%LTSOLVE Utility routine for "preconditioned" iterative methods.
%
% x = ltsolve(L,y,W,T)
%
% Computes the vector
%   x = (L_p)'*y
% where L_p is the A-weighted generalized inverse of L.
%
% Here, L is a p-by-n matrix, W holds a basis for the null space of L,
% and T is a utility matrix which should be computed by routine pinit.
%
% Alternatively, L is square, and W and T are not needed.
%
% If W and T are not specified, then instead the routine computes
%   x = inv(L(:,1:p))'*y(1:p) .
%
% Notice that x and y may be matrices, in which case x(:,i)
% corresponds to y(:,i).

% Reference: P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, 07/29/97.

% Initialization.
[p,n] = size(L); nu = n-p;

% Special treatment of square L.
if (nu==0), x = (L')\y; return; end

% Perform the first stage, if necessary.
if (nargin > 2), y = y(1:p) - T(:,1:p)'*(W'*y); end

% Always perform the second stage.
x = y(1:p)'/L(:,1:p); x = x';