function [X,rho,eta] = pmr2(A,L,N,b,k,reorth)
%PMR2 Preconditioned MR-II algorithm for symmetric indefinite problems
%
% [X,rho,eta] = pmr2(A,L,N,b,k,reorth)
%
% PMR2 applies smoothing-norm preconditioning to the MR-II method, which
% is a variant the MINRES algorithm for symmetric indefinite linear
% systems A x = b, with starting vector A*b (instead of b as in MINRES).
% This function returns all k iterates, stored as the columns of the
% matrix X.  The solution norm and residual norm are returned in eta and
% rho, respectively.
%
% The preconditioner uses two matrices: the matrix L that defines the
% smoothing norm, and the matrix N whose columns span the null space
% of L.  It is assumed that L is p-times-n with p < n.
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% Reference: P. C. Hansen and T. K. Jensen, "Smoothing-norm preconditioning
% for regularizing minimum-residual methods", SIAM J. Matrix Anal. Appl.
% 29 (2006), 1-14.

% Per Christian Hansen, IMM, September 21, 2007.

% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==5), reorth = 0; end
[m,n] = size(A);
p = size(L,1);
if (m ~= n || norm(A-A','fro')), error('The matrix must be symmetric'), end

% Allocate space.
X = zeros(n,k);
if reorth
  W = zeros(p,k);
  if (k>=n), error('No. of iterations must satisfy k < n'), end
end
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
end

% Initialization for working with pseudoinverses of L.
[Q0,R0] = qr(A*N,0);  % Compate QR factorization of A*N.
T0 = N'*Q0;
[QL,RL] = qr(L',0);   % Compact QR factgorization of L'.
TN = pinit(N,A);      % Prepare for A-weighted pseudoinverse computations.
bb = RL\( QL'*( b - Q0*(T0\(N'*b)) ) );

% Prepare for interation.
x0 = N*( R0\(Q0'*b) );
xi = zeros(p,1); r = bb;
vold = 0;
v1 = A*( QL*(RL'\r) );
v2 = Q0*( T0\(N'*v1) );
v  = RL\( QL'*(v1-v2) );
wold = 0;
v1 = A*( QL*(RL'\v) );
v2 = Q0*( T0\(N'*v1) );
w  = RL\( QL'*(v1-v2) );
beta = norm(w);
v = v./beta; w = w./beta;
if reorth, W(:,1) = w; end

% Perform k iterations.
for i=1:k
    
  rrho = r'*w;
  xi = xi + rrho*v;
  r  = r  - rrho*w;
  v1 = A*( QL*(RL'\w) );
  v2 = Q0*( T0\(N'*v1) );
  Aw = RL\( QL'*(v1-v2) );
  alpha = w'*Aw;
  vnew  =  w - alpha*v - beta*vold;
  wnew  = Aw - alpha*w - beta*wold;
  vold = v; wold = w; v = vnew; w = wnew;
  if reorth
      for j=1:i, w = w - (W(:,j)'*w)*W(:,j); end
  end;
  beta = norm(w);
  v = v./beta;	w = w./beta;
  if reorth, W(:,i+1) = w; end;

  X(:,i) = lsolve(L,xi,N,TN) + x0;
  if (nargout>1), rho(i) = norm(A*X(:,i)-b); end
  if (nargout>2), eta(i) = norm(xi); end
  
end