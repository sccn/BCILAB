function [A,U,V] = regutm(m,n,s)
%REGUTM Test matrix for regularization methods.
%
% [A,U,V] = regutm(m,n,s)
%
% Generates a random m-times-n matrix A such that A*A' and A'*A
% are oscillating.  Hence, in the SVD of A,
%    A = U*diag(s)*V',
% the number of sign changes in U(:,i) and V(:,i) is exactly i-1.
%
% The third argument s specifies the singular values of A.  If not
% present, then s = logspace(0,round(log10(eps)),n).

% Reference: P. C. Hansen, "Test matrices for regularization methods",
% SIAM J. Sci. Comput. 16 (1995), 506--512.

% Per Christian Hansen, IMM, 07/30/97.

% Initialization.
if (nargin==1), n = m; end
if (nargin<3), s = logspace(0,round(log10(eps)),min(m,n)); end

% Special treatment of the case m < n.
if (m < n), [A,V,U] = regutm(n,m,s); A = A'; return, end

% Generate random bidiagonal matrix with nonnegative elements.
if (n < 100), mu = .222*n + .0278*n^2; else mu = 3*n; end
B = abs(diag(randn(n,1)+mu) + diag(randn(n-1,1)+mu,1));

% Compute the SVD of B.
[U,dummy,V] = svd(B); clear dummy

% Repeat if m > n.
if (m > n)
  clear U
  B = abs(diag(randn(m,1)+mu) + diag(randn(m-1,1)+mu,1));
  [U,dummy] = svd(B); clear dummy, U = U(:,1:n);
end

% Compute A.
A = U*diag(s)*V';