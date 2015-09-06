function [X,rho,eta,F] = pcgls(A,L,W,b,k,reorth,sm)
%PCGLS "Precond." conjugate gradients appl. implicitly to normal equations.
% [X,rho,eta,F] = pcgls(A,L,W,b,k,reorth,sm)
%
% Performs k steps of the `preconditioned' conjugate gradient
% algorithm applied implicitly to the normal equations
%    (A*L_p)'*(A*L_p)*x = (A*L_p)'*b ,
% where L_p is the A-weighted generalized inverse of L.  Notice that the
% matrix W holding a basis for the null space of L must also be specified.
%
% The routine returns all k solutions, stored as columns of the matrix X.
% The solution seminorm and residual norm are returned in eta and rho,
% respectively.
%
% If the generalized singular values sm of (A,L) are also provided,
% pcgls computes the filter factors associated with each step and
% stores them columnwise in the matrix F.
%
% Reorthogonalization of the normal equation residual vectors
% A'*(A*X(:,i)-b) is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% References: A. Bjorck, "Numerical Methods for Least Squares Problems",
% SIAM, Philadelphia, 1996.
% P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM and Martin Hanke, Institut fuer
% Praktische Mathematik, Universitaet Karlsruhe, 07/02/97.

% The fudge threshold is used to prevent filter factors from exploding.
fudge_thr = 1e-4;

% Initialization
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==5), reorth = 0; end
if (nargout==4 & nargin<7), error('Too few input arguments'), end
if (reorth<0 | reorth>1), error('Illegal reorth'), end
[m,n] = size(A); p = size(L,1); X = zeros(n,k);
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
end
if (nargin==7)
  F = zeros(p,k); Fd = zeros(p,1); gamma = (sm(:,1)./sm(:,2)).^2;
end

% Prepare for computations with L_p.
[NAA,x_0] = pinit(W,A,b);

% Prepare for CG iteartion.
x  = x_0;
r  = b - A*x_0; s = A'*r;
q1 = ltsolve(L,s);
q  = lsolve(L,q1,W,NAA);
z  = q;
dq = s'*q;
if (nargout>2), z1 = q1; x1 = zeros(p,1); end
if (reorth==1)
  Q1n = zeros(p,k);
  Q1n = q1/norm(q1);
end

% Iterate.
for j=1:k

  % Update x and r vectors; compute q1.
  Az  = A*z; alpha = dq/(Az'*Az);
  x   = x + alpha*z;
  r   = r - alpha*Az; s = A'*r;
  q1  = ltsolve(L,s);

  % Reorthogonalize q1 to previous q1-vectors, if required.
  if (reorth==1)
    for i=1:j, q1 = q1 - (Q1n(:,i)'*q1)*Q1n(:,i); end
    Q1n = [Q1n,q1/norm(q1)];
  end

  % Update z vector.
  q   = lsolve(L,q1,W,NAA);
  dq2 = s'*q; beta = dq2/dq;
  dq  = dq2;
  z   = q + beta*z;
  X(:,j) = x;
  if (nargout>1), rho(j) = norm(r); end
  if (nargout>2)
    x1 = x1 + alpha*z1; z1 = q1 + beta*z1; eta(j) = norm(x1);
  end

  % Compute filter factors, if required.
  if (nargin==7)
    if (j==1)
      F(:,1) = alpha*gamma;
      Fd = gamma - gamma.*F(:,1) + beta*gamma;
    else
      F(:,j) = F(:,j-1) + alpha*Fd;
      Fd = gamma - gamma.*F(:,j) + beta*Fd;
    end
    if (j > 2)
      f = find(abs(F(:,j-1)-1) < fudge_thr & abs(F(:,j-2)-1) < fudge_thr);
      if ~isempty(f), F(f,j) = ones(length(f),1); end
    end
  end

end