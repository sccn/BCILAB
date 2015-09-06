function x = lsolve(L,y,W,T)
%LSOLVE Utility routine for "preconditioned" iterative methods.
%
% x = lsolve(L,y,W,T)
%
% Computes the vector
%    x = L_p*y
% where L_p is the A-weighted generalized inverse of L.
%
% Here, L is a p-by-n matrix, W holds a basis for the null space of L,
% and T is a utility matrix which should be computed by routine pinit.
%
% Alternatively, L is square, and W and T are not needed.
%
% Notice that x and y may be matrices, in which case
%    x(:,i) = L_p*y(:,i) .

% Reference: P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, 07/29/97.

% Initialization.
[p,n] = size(L); nu = n-p; ly = size(y,2);

% Special treatment of square L.
if (nu==0), x = L\y; return; end

% The general case.
x = L(:,1:p)\y;
x = [x;zeros(nu,ly)] - W*(T(:,1:p)*x);