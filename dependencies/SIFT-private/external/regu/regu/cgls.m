function [X,rho,eta,F] = cgls(A,b,k,reorth,s)
%CGLS Conjugate gradient algorithm applied implicitly to the normal equations.
%
% [X,rho,eta,F] = cgls(A,b,k,reorth,s)
%
% Performs k steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The corresponding solution and residual norms
% are returned in the vectors eta and rho, respectively.
%
% If the singular values s are also provided, cgls computes the
% filter factors associated with each step and stores them
% columnwise in the matrix F.
%
% Reorthogonalization of the normal equation residual vectors
% A'*(A*X(:,i)-b) is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% References: A. Bjorck, "Numerical Methods for Least Squares Problems",
% SIAM, Philadelphia, 1996.
% C. R. Vogel, "Solving ill-conditioned linear systems using the
% conjugate gradient method", Report, Dept. of Mathematical
% Sciences, Montana State University, 1987.

% Per Christian Hansen, IMM, July 23, 2007.

% The fudge threshold is used to prevent filter factors from exploding.
fudge_thr = 1e-4;

% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==3), reorth = 0; end
if (nargout==4 & nargin<5), error('Too few input arguments'), end
if (reorth<0 | reorth>1), error('Illegal reorth'), end
[m,n] = size(A); X = zeros(n,k);
if (reorth==1), ATr = zeros(n,k+1); end
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
end
if (nargin==5)
  F = zeros(n,k); Fd = zeros(n,1); s2 = s.^2;
end

% Prepare for CG iteration.
x = zeros(n,1);
d = A'*b;
r = b;
normr2 = d'*d;
if (reorth==1), ATr(:,1) = d/norm(d); end

% Iterate.
for j=1:k

  % Update x and r vectors.
  Ad = A*d; alpha = normr2/(Ad'*Ad);
  x  = x + alpha*d;
  r  = r - alpha*Ad;
  s  = A'*r;

  % Reorthogonalize s to previous s-vectors, if required.
  if (reorth==1)
    for i=1:j, s = s - (ATr(:,i)'*s)*ATr(:,i); end
    ATr(:,j+1) = s/norm(s);
  end

  % Update d vector.
  normr2_new = s'*s;
  beta = normr2_new/normr2;
  normr2 = normr2_new;
  d = s + beta*d;
  X(:,j) = x;

  % Compute norms, if required.
  if (nargout>1), rho(j) = norm(r); end
  if (nargout>2), eta(j) = norm(x); end

  % Compute filter factors, if required.
  if (nargin==5)
    if (j==1)
      F(:,1) = alpha*s2;
      Fd = s2 - s2.*F(:,1) + beta*s2;
    else
      F(:,j) = F(:,j-1) + alpha*Fd;
      Fd = s2 - s2.*F(:,j) + beta*Fd;
    end
    if (j > 2)
      f = find(abs(F(:,j-1)-1) < fudge_thr & abs(F(:,j-2)-1) < fudge_thr);
      if ~isempty(f), F(f,j) = ones(length(f),1); end
    end
  end

end