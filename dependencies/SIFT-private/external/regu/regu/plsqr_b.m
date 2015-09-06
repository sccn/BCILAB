function [X,rho,eta,F] = plsqr_b(A,L,W,b,k,reorth,sm)
%PLSQR_B "Precond." version of the LSQR Lanczos bidiagonalization algorithm.
%
% [X,rho,eta,F] = plsqr_b(A,L,W,b,k,reorth,sm)
%
% Performs k steps of the `preconditioned' LSQR Lanczos
% bidiagonalization algorithm applied to the system
%    min || (A*L_p) x - b || ,
% where L_p is the A-weighted generalized inverse of L.  Notice
% that the matrix W holding a basis for the null space of L must
% also be specified.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The solution seminorm and the residual norm are
% returned in eta and rho, respectively.
%
% If the generalized singular values sm of (A,L) are also provided,
% then glsqr computes the filter factors associated with each step
% and stores them columnwise in the matrix F.
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS

% References: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
% sparse linear equations and sparse least squares", ACM Trans.
% Math. Software 8 (1982), 43-71.
% P. C. Hansen, "Rank-Deficient and Discrete Ill-Posed Problems.
% Numerical Aspects of Linear Inversion", SIAM, Philadelphia, 1997.

% Per Christian Hansen, IMM, Sept. 13, 2001.

% The fudge threshold is used to prevent filter factors from exploding.
fudge_thr = 1e-4;

% Initialization
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==5), reorth = 0; end
if (nargout==4 & nargin<7), error('Too few input arguments'), end
[m,n] = size(A); X = zeros(n,k); [pp,n1] = size(L);
if (n1 ~= n | m < n | n < pp)
  error('Incorrect dimensions of A and L')
end
if (reorth==0)
  UV = 0;
elseif (reorth==1)
  if (k>=n), error('No. of iterations must satisfy k < n'), end
  U = zeros(m,k); V = zeros(pp,k); UV = 1;
else
  error('Illegal reorth')
end
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
  c2 = -1; s2 = 0; xnorm = 0; z = 0;
end
if (nargin==7)
  ls = size(sm,1);
  F = zeros(ls,k); Fv = zeros(ls,1); Fw = Fv;
  s = (sm(:,1)./sm(:,2)).^2;
end

% Prepare for computations with L_p.
[NAA,x_0] = pinit(W,A,b);

% By subtracting the component A*x_0 from b we insure that
% the corrent residual norms are computed.
b = b - A*x_0;

% Prepare for LSQR iteration.
v = zeros(pp,1); x = v; beta = norm(b);
if (beta==0), error('Right-hand side must be nonzero'), end
u = b/beta; if (UV), U(:,1) = u; end
r = ltsolve(L,A'*u,W,NAA); alpha = norm(r);
v = r/alpha; if (UV), V(:,1) = v; end
phi_bar = beta; rho_bar = alpha; w = v;
if (nargin==7), Fv = s/(alpha*beta); Fw = Fv; end

% Perform Lanczos bidiagonalization with/without reorthogonalization.
for i=2:k+1

  alpha_old = alpha; beta_old = beta;

  % Compute (A*L_p)*v - alpha*u.
  p = A*lsolve(L,v,W,NAA) - alpha*u;
  if (reorth==0)
    beta = norm(p); u = p/beta;
  else
    for j=1:i-1, p = p - (U(:,j)'*p)*U(:,j); end
    beta = norm(p); u = p/beta;
  end

  % Compute L_p'*A'*u - beta*v.
  r = ltsolve(L,A'*u,W,NAA) - beta*v;
  if (reorth==0)
    alpha = norm(r); v = r/alpha;
  else
    for j=1:i-1, r = r - (V(:,j)'*r)*V(:,j); end
    alpha = norm(r); v = r/alpha;
  end

  % Store U and V if necessary.
  if (UV), U(:,i) = u; V(:,i) = v; end

  % Construct and apply orthogonal transformation.
  rrho = norm([rho_bar,beta]); c1 = rho_bar/rrho;
  s1 = beta/rrho; theta = s1*alpha; rho_bar = -c1*alpha;
  phi = c1*phi_bar; phi_bar = s1*phi_bar;

  % Compute solution norm and residual norm if necessary;
  if (nargout > 1)
    delta = s2*rrho; gamma_bar = -c2*rrho; rhs = phi - delta*z;
    z_bar = rhs/gamma_bar; eta(i-1) = norm([xnorm,z_bar]);
    gamma = norm([gamma_bar,theta]);
    c2 = gamma_bar/gamma; s2 = theta/gamma;
    z = rhs/gamma; xnorm = norm([xnorm,z]);
    rho(i-1) = abs(phi_bar);
  end

  % If required, compute the filter factors.
  if (nargin==7)

    if (i==2)
      Fv_old = Fv;
      Fv  = Fv.*(s - beta^2 - alpha_old^2)/(alpha*beta);
      F(:,i-1) = (phi/rrho)*Fw;
    else
      tmp = Fv;
      Fv = (Fv.*(s - beta^2 - alpha_old^2) - ...
                 Fv_old*alpha_old*beta_old)/(alpha*beta);
      Fv_old = tmp;
      F(:,i-1) = F(:,i-2) + (phi/rrho)*Fw;
    end
    if (i > 3)
      f = find(abs(F(:,i-2)-1) < fudge_thr & abs(F(:,i-3)-1) < fudge_thr);
      if ~isempty(f), F(f,i-1) = ones(length(f),1); end
    end
    Fw = Fv - (theta/rrho)*Fw;

  end

  % Update the solution.
  x = x + (phi/rrho)*w; w = v - (theta/rrho)*w;
  X(:,i-1) = lsolve(L,x,W,NAA) + x_0;

end