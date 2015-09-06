function [U,B_k,V] = lanc_b(A,p,k,reorth)
%LANC_B Lanczos bidiagonalization.
%
% B_k = lanc_b(A,p,k,reorth)
% [U,B_k,V] = lanc_b(A,p,k,reorth)
%
% Performs k steps of the Lanczos bidiagonalization process with
% starting vector p, producing a lower bidiagonal matrix
%           [b_11               ]
%           [b_21 b_22          ]
%     B_k = [     b_32 .        ]
%           [          . b_kk   ]
%           [            b_k+1,k]
% such that
%     A*V = U*B_k ,
% where U and V consist of the left and right Lanczos vectors.
%
% Reorthogonalization is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization,
%    reorth = 1 : reorthogonalization by means of MGS,
%    reorth = 2 : Householder-reorthogonalization.
% No reorthogonalization is assumed if reorth is not specified.

% Reference: G. H. Golub & C. F. Van Loan, "Matrix Computations",
% 3. Ed., Johns Hopkins, 1996.  Section 9.3.4.
% Referred to as "bidiag1" by Paige and Saunders.

% Per Christian Hansen, IMM, April 8, 2001.

% Initialization.
if (k<1), error('Number of steps k must be positive'), end
if (nargin < 4), reorth = 0; end
if (reorth < 0 | reorth > 2), error('Illegal reorth'), end
if (nargout==2), error('Not enough output arguments'), end
[m,n] = size(A);
B_k = sparse(k+1,k);
if (nargout>1 | reorth==1)
  U = zeros(m,k); V = zeros(n,k); UV = 1;
else
  UV = 0;
end
if (reorth==2)
  if (k>=n), error('No. of iterations must satisfy k < n'), end
  HHU = zeros(m,k); HHV = zeros(n,k);
  HHalpha = zeros(1,k); HHbeta = HHalpha;
end

% Prepare for Lanczos iteration.
v = zeros(n,1);
beta = norm(p);
if (beta==0), error('Starting vector must be nonzero'), end
if (reorth==2)
  [beta,HHbeta(1),HHU(:,1)] = gen_hh(p);
end
u = p/beta;
if (UV), U(:,1) = u; end

% Perform Lanczos bidiagonalization with/without reorthogonalization.
for i=1:k

  r = A'*u - beta*v;
  if (reorth==0)
    alpha = norm(r); v = r/alpha;
  elseif (reorth==1)
    for j=1:i-1, r = r - (V(:,j)'*r)*V(:,j); end
    alpha = norm(r); v = r/alpha;
  else
    for j=1:i-1
      r(j:n) = app_hh(r(j:n),HHalpha(j),HHV(j:n,j));
    end
    [alpha,HHalpha(i),HHV(i:n,i)] = gen_hh(r(i:n));
    v = zeros(n,1); v(i) = 1;
    for j=i:-1:1
      v(j:n) = app_hh(v(j:n),HHalpha(j),HHV(j:n,j));
    end
  end
  B_k(i,i) = alpha; if (UV), V(:,i) = v; end

  p = A*v - alpha*u;
  if (reorth==0)
    beta = norm(p); u = p/beta;
  elseif (reorth==1)
    for j=1:i, p = p - (U(:,j)'*p)*U(:,j); end
    beta = norm(p); u = p/beta;
  else
    for j=1:i
      p(j:m) = app_hh(p(j:m),HHbeta(j),HHU(j:m,j));
    end
    [beta,HHbeta(i+1),HHU(i+1:m,i+1)] = gen_hh(p(i+1:m));
    u = zeros(m,1); u(i+1) = 1;
    for j=i+1:-1:1
      u(j:m) = app_hh(u(j:m),HHbeta(j),HHU(j:m,j));
    end
  end
  B_k(i+1,i) = beta; if (UV), U(:,i+1) = u; end

end

if (nargout==1), U = B_k; end