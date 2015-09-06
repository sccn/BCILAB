function [x_k,rho,eta] = tsvd(U,s,V,b,k)
%TSVD Truncated SVD regularization.
%
% [x_k,rho,eta] = tsvd(U,s,V,b,k)
%
% Computes the truncated SVD solution
%    x_k = V(:,1:k)*inv(diag(s(1:k)))*U(:,1:k)'*b .
% If k is a vector, then x_k is a matrix such that
%    x_k = [ x_k(1), x_k(2), ... ] .
%
% The solution and residual norms are returned in eta and rho.

% Per Christian Hansen, IMM, 12/21/97.

% Initialization.
[n,p] = size(V); lk = length(k);
if (min(k)<0 | max(k)>p)
  error('Illegal truncation parameter k')
end
x_k = zeros(n,lk);
eta = zeros(lk,1); rho = zeros(lk,1);
beta = U(:,1:p)'*b;
xi = beta./s;

% Treat each k separately.
for j=1:lk
  i = k(j);
  if (i>0)
    x_k(:,j) = V(:,1:i)*xi(1:i);
    eta(j) = norm(xi(1:i));
    rho(j) = norm(beta(i+1:p));
  end
end

if (nargout > 1 & size(U,1) > p)
  rho = sqrt(rho.^2 + norm(b - U(:,1:p)*beta)^2);
end