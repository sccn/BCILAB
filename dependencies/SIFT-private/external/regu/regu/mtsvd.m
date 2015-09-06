function [x_k,rho,eta] = mtsvd(U,s,V,b,k,L)
%MTSVD Modified truncated SVD regularization.
%
% [x_k,rho,eta] = mtsvd(U,s,V,b,k,L)
%
% Computes the modified TSVD solution:
%    x_k = V*[ xi_k ] .
%            [ xi_0 ]
% Here, xi_k defines the usual TSVD solution
%    xi_k = inv(diag(s(1:k)))*U(:,1:k)'*b ,
% and xi_0 is chosen so as to minimize the seminorm || L x_k ||.
% This leads to choosing xi_0 as follows:
%    xi_0 = -pinv(L*V(:,k+1:n))*L*V(:,1:k)*xi_k .
%
% The truncation parameter must satisfy k > n-p.
%
% If k is a vector, then x_k is a matrix such that
%     x_k = [ x_k(1), x_k(2), ... ] .
%
% The solution and residual norms are returned in eta and rho.

% Reference: P. C. Hansen, T. Sekii & H. Shibahashi, "The modified
% truncated-SVD method for regularization in general form", SIAM J.
% Sci. Stat. Comput. 13 (1992), 1142-1150.

% Per Christian Hansen, IMM, 12/22/95.

% Initialization.
m = size(U,1); [p,n] = size(L);
lk = length(k); kmin = min(k);
if (kmin<n-p+1 | max(k)>n)
  error('Illegal truncation parameter k')
end
x_k = zeros(n,lk);
beta = U(:,1:n)'*b; xi = beta./s;
eta = zeros(lk,1); rho =zeros(lk,1);

% Compute large enough QR factorization.
[Q,R] = qr(L*V(:,n:-1:kmin+1),0);

% Treat each k separately.
for j=1:lk
  kj = k(j); xtsvd = V(:,1:kj)*xi(1:kj);
  if (kj==n)
    x_k(:,j) = xtsvd;
  else
    z = R(1:n-kj,1:n-kj)\(Q(:,1:n-kj)'*(L*xtsvd));
    z = z(n-kj:-1:1);
    x_k(:,j) = xtsvd - V(:,kj+1:n)*z;
  end
  eta(j) = norm(x_k(:,j));
  rho(j) = norm(beta(kj+1:n) + s(kj+1:n).*z);
end

if (nargout > 1 & m > n)
  rho = sqrt(rho.^2 + norm(b - U(:,1:n)*beta)^2);
end