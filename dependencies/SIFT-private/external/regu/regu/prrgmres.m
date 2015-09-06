function [X,rho,eta] = prrgmres(A,L,N,b,k)
%PRRGMRES Preconditioned RRGMRES algorithm for square inconsistent systems
%
% [X,rho,eta] = rrgmres(A,L,N,b,k)
%
% PRRGMRES applies smoothing-norm preconditioning to the RRGMRES method,
% which is a variant of RRGMRES for square linear systems A x = b, with  
% starting vector A*b (instead of b as in GMRES).  This function returns
% all k iterates, stored as columns of the matrix X.  The solution norm
% and residual norm are returned in eta and rho, resp.
%
% The preconditioner uses two matrices: the matrix L that defines the
% smoothing norm, and the matrix N whose columns span the null space
% of L.  It is assumed that L is p-times-n with p < n.
%
% For symmetric matrices, use the function pmr2 instead.

% Reference: P. C. Hansen and T. K. Jensen, "Smoothing-norm preconditioning
% for regularizing minimum-residual methods", SIAM J. Matrix Anal. Appl.
% 29 (2006), 1-14.

% Per Christian Hansen, IMM, September 21, 2007.

% Check input arguments.
if (k < 1), error('Number of steps k must be positive'), end
[m,n] = size(A);
p = size(L,1);
if (m ~= n), error('A must be square'), end
if nargin < 5, error('Too few input arguments'), end

% Allocate space.
X = zeros(n,k); % Matrix of solutions.
V = zeros(p,k); % Orthonormal vectors spanning Krylov subspace
h = zeros(k,1); % New column of Hessenberg matrix
Q = zeros(k+1); % H = Q*T, Q orthogonal
T = zeros(k);   %          T upper triangular
W = zeros(p,k); % W = V*inv(T)
if (nargout>1), rho = zeros(k,1); end
if (nargout>2), eta = zeros(k,1); end

% Initialization for working with pseudoinverses of L.
[Q0,R0] = qr(A*N,0);  % Compate QR factorization of A*N.
T0 = N'*Q0;
[QL,RL] = qr(L',0);   % Compact QR factgorization of L'.
TN = pinit(N,A);      % Prepare for A-weighted pseudoinverse computations.
bb = RL\( QL'*( b - Q0*(T0\(N'*b)) ) );

% Initialize variables.
v1     = A*( QL*(RL'\bb) );
v2     = Q0*( T0\(N'*v1) );
r      = RL\( QL'*(v1-v2) );
alpha  = norm(r);
V(:,1) = r/alpha;  % Initial vector of Krylov subspace is A*b.
Q(1,1) = 1;
x0     = N*( R0\(Q0'*b) );
xi     = zeros(p,1);
beta   = V(:,1)'*bb;

% Begin iterations.
for i=1:k

  v1 = A*( QL*(RL'\V(:,i)) );
  v2 = Q0*( T0\(N'*v1) );
  r  = RL\( QL'*(v1-v2) );

  % Modified Gram-Schmidt on the new vector.
  for j=1:k
    h(j) = V(:,j)'*r;
    r = r - V(:,j)*h(j);
  end
  alpha = norm(r);

  % Store new Arnoldi vector and update projected rhs.
  V(:,i+1) = r/alpha;
  beta = [beta;V(:,i+1)'*bb];

  % Apply previous rotations to h.
  T(1:i,i) = Q(1:i,1:i)'*h(1:i);

  % Compute Givens rotation parameters.
  rc = T(i,i);
  if alpha == 0
    c = 1; s = 0;
  elseif abs(alpha) > abs(rc)
    tau = -rc/alpha;
    s = 1 / sqrt(1 + abs(tau)^2);
    c = s*tau;
  else
    tau = -alpha/rc;
    c = 1 / sqrt(1 + abs(tau)^2);
    s = c*tau;
  end

  % Apply givens rotations.
  T(i,i) = c'*rc - s'*alpha;
  Q(1:i,[i,i+1]) = Q(1:i,i)*[c s];
  Q(i+1,[i,i+1]) = [-s c];
    
  if abs(T(i,i)) <= eps
    disp('Hession matrix is (numerically) singular')
  end

  % Update W = V*inv(T);
  W(:,i) = (V(:,i) - W(:,1:i-1)*T(1:i-1,i))/T(i,i);

  % Update solution.
  xi = xi + (Q(1:i+1,i)'*beta)*W(:,i);

  % Update output variables.
  X(:,i) = lsolve(L,xi,N,TN) + x0;
  if nargout>1, rho(i) = norm(A*X(:,i)-b); end
  if nargout>2, eta(i) = norm(xi); end

end