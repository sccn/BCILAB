function [X,rho,eta] = mr2(A,b,k,reorth)
%MR2 Solution of symmetric indefinite problems by the MR-II algorithm
%
% [X,rho,eta] = mr2(A,b,k,reorth)
%
% MR-II is a variant of the MINRES algorithm for symmetric indefinite linear
% systems A x = b, with starting vector A*b (instead of b as in MINRES).
% This function returns all k iterates, stored as the columns of the
% matrix X.  The solution norm and residual norm are returned in eta and
% rho, respectively.
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% Reference: M. Hanke, "Conjugate Gradient Methods for Ill-Posed Problems",
% Longman Scientific and Technical, Essex, 1995.

% Per Christian Hansen, IMM, September 1, 2007.
% Based on the function mr2 from Restore Tools by James G. Nagy.

% Initialization.
if (k < 1), error('Number of steps k must be positive'), end
if (nargin==3), reorth = 0; end
[m,n] = size(A);
if (m ~= n || norm(A-A','fro')), error('The matrix must be symmetric'), end

% Allocate space.
X = zeros(n,k);
if reorth
  W = zeros(n,k);
  if (k>=n), error('No. of iterations must satisfy k < n'), end
end
if (nargout > 1)
  eta = zeros(k,1); rho = eta;
end

% Prepare for interation.
x = zeros(n,1); r = b;
vold = 0; v = A*r;
wold = 0; w = A*v;
beta = norm(w);
v = v./beta; w = w./beta;
if reorth, W(:,1) = w; end

% Perform k iterations.
for i=1:k
    
  rrho = r'*w;
  x = x + rrho*v;
  r = r - rrho*w;
  Aw = A*w;
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

  X(:,i) = x;
  if (nargout>1), rho(i) = norm(r); end
  if (nargout>2), eta(i) = norm(x); end
  
end