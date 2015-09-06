function [X,rho,eta,F] = nu(A,b,k,nu,s)
%NU Brakhage's nu-method.
%
% [X,rho,eta,F] = nu(A,b,k,nu,s)
%
% Performs k steps of Brakhage's nu-method for the problem
%    min || A x - b || .
% The routine returns all k solutions, stored as columns of
% the matrix X.  The solution norm and residual norm are returned
% in eta and rho, respectively.
%
% If nu is not specified, nu = .5 is the default value, which gives
% the Chebychev method of Nemirovskii and Polyak.
%
% If the singular values s are also provided, nu computes the
% filter factors associated with each step and stores them
% columnwise in the matrix F.

% Reference: H. Brakhage, "On ill-posed problems and the method of
% conjugate gradients"; in H. W. Engl & G. W. Groetsch, "Inverse and
% Ill-Posed Problems", Academic Press, 1987.

% Martin Hanke, Institut fuer Praktische Mathematik, Universitaet
% Karlsruhe and Per Christian Hansen, IMM, 03/21/92.

% Set parameter.
l_steps = 3;      % Number of Lanczos steps for est. of || A ||.
fudge   = 0.99;   % Scale A and b by fudge/|| A*L_p ||.
fudge_thr = 1e-4; % Used to prevent filter factors from exploding.

% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==3), nu = .5; end
[m,n] = size(A); X = zeros(n,k);
if (nargout > 1)
  rho = zeros(k,1); eta = rho;
end;
if (nargin==5)
  F = zeros(n,k); Fd = zeros(n,1); s = s.^2;
end
V = zeros(n,l_steps); B = zeros(l_steps+1,l_steps);
v = zeros(n,1); eta = zeros(l_steps+1,1);

% Compute a rough estimate of the norm of A by means of a few
% steps of Lanczos bidiagonalization, and scale A and b such
% that || A || is slightly less than one.
beta = norm(b); u = b/beta;
for i=1:l_steps
  r = A'*u - beta*v;
  alpha = norm(r); v = r/alpha;
  B(i,i) = alpha; V(:,i) = v;
  p = A*v - alpha*u;
  beta = norm(p); u = p/beta;
  B(i+1,i) = beta;
end
scale = fudge/norm(B); A = scale*A; b = scale*b;
if (nargin==5), s = scale^2*s; end

% Prepare for iteration.
x = zeros(n,1);
d = A'*b;
r = d;
if (nargout>1), z = b; end

% Iterate.
for j=0:k-1

  % Updates.
  alpha = 4*(j+nu)*(j+nu+0.5)/(j+2*nu)/(j+2*nu+0.5);
  beta  = (j+nu)*(j+1)*(j+0.5)/(j+2*nu)/(j+2*nu+0.5)/(j+nu+1);
  Ad = A*d; AAd = A'*Ad;
  x  = x + alpha*d;
  r  = r - alpha*AAd;
  d  = r + beta*d;
  X(:,j+1) = x;
  if (nargout>1)
    z = z - alpha*Ad; rho(j+1) = norm(z)/scale;
  end;
  if (nargout>2), eta(j+1) = norm(x); end;

  % Filter factors.
  if (nargin==5)
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