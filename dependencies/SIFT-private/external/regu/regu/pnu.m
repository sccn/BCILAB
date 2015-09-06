function [X,rho,eta,F] = pnu(A,L,W,b,k,nu,sm)
%PNU "Preconditioned" version of Brakhage's nu-method.
%
% [X,rho,eta,F] = pnu(A,L,W,b,k,nu,sm)
%
% Performs k steps of a `preconditioned' version of Brakhage's
% nu-method for the problem
%    min || (A*L_p) x - b || ,
% where L_p is the A-weighted generalized inverse of L.  Notice
% that the matrix W holding a basis for the null space of L must
% also be specified.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The solution seminorm and residual norm are returned
% in eta and rho, respectively.
%
% If nu is not specified, nu = .5 is the default value, which gives
% the Chebychev method of Nemirovskii and Polyak.
%
% If the generalized singular values sm of (A,L) are also provided,
% then pnu computes the filter factors associated with each step and
% stores them columnwise in the matrix F.

% Reference: H. Brakhage, "On ill-posed problems and the method of
% conjugate gradients"; in H. W. Engl & G. W. Groetsch, "Inverse and
% Ill-Posed Problems", Academic Press, 1987.

% Martin Hanke, Institut fuer Praktische Mathematik, Universitaet
% Karlsruhe and Per Christian Hansen, IMM, 06/25/92.

% Set parameters.
l_steps = 3;      % Number of Lanczos steps for est. of || A*L_p ||.
fudge   = 0.99;   % Scale A and b by fudge/|| A*L_p ||.
fudge_thr = 1e-4; % Used to prevent filter factors from exploding.
 
% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==5), nu = .5; end
[m,n] = size(A); p = size(L,1); X = zeros(n,k);
if (nargout > 1)
  rho = zeros(k,1); eta = rho;
end;
if (nargin==7)
  F = zeros(n,k); Fd = zeros(n,1); s = (sm(:,1)./sm(:,2)).^2;
end
V = zeros(p,l_steps); B = zeros(l_steps+1,l_steps);
v = zeros(p,1);

% Prepare for computations with L_p.
[NAA,x_0] = pinit(W,A,b); x1 = x_0;

% Compute a rough estimate of || A*L_p || by means of a few steps of
% Lanczos bidiagonalization, and scale A and b such that || A*L_p || is
% slightly less than one.
b_0 = b - A*x_0; beta = norm(b_0); u = b_0/beta;
for i=1:l_steps
  r = ltsolve(L,A'*u,W,NAA) - beta*v;
  alpha = norm(r); v = r/alpha;
  B(i,i) = alpha; V(:,i) = v;
  p = A*lsolve(L,v,W,NAA) - alpha*u;
  beta = norm(p); u = p/beta;
  B(i+1,i) = beta;
end
scale = fudge/norm(B); A = scale*A;
if (nargin==7), s = scale^2*s; end

% Prepare for iteration.
x  = x_0;
z  = -scale*b_0;
r  = A'*z;
d1 = ltsolve(L,r);
d  = lsolve(L,d1,W,NAA);
if (nargout>2), x1 = L*x_0; end

% Iterate.
for j=0:k-1

  % Updates.
  alpha = 4*(j+nu)*(j+nu+0.5)/(j+2*nu)/(j+2*nu+0.5);
  beta  = -(j+nu)*(j+1)*(j+0.5)/(j+2*nu)/(j+2*nu+0.5)/(j+nu+1);
  Ad  = A*d; AAd = A'*Ad;
  x   = x - alpha*d;
  r   = r - alpha*AAd;
  rr1 = ltsolve(L,r);
  rr  = lsolve(L,rr1,W,NAA);
  d   = rr - beta*d;
  X(:,j+1) = x;
  if (nargout>1 )
    z = z - alpha*Ad; rho(j+1) = norm(z)/scale;
  end
  if (nargout>2)
    x1 = x1 - alpha*d1; d1 = rr1 - beta*d1;
    eta(j+1) = norm(x1);
  end

  % Filter factors.
  if (nargin==7)
    if (j==0)
      F(:,1) = alpha*s;
      Fd = s - s.*F(:,1) + beta*s;
    else
      F(:,j+1) = F(:,j) + alpha*Fd;
      Fd = s - s.*F(:,j+1) + beta*Fd;
    end
    if (j > 1)
      f = find(abs(F(:,j)-1) < fudge_thr & abs(F(:,j-1)-1) < fudge_thr);
      if ~isempty(f), F(f,j+1) = ones(length(f),1); end
    end
  end

end