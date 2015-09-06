function [U,B,V] = bidiag(A)
%BIDIAG Bidiagonalization of an m-times-n matrix with m >= n.
%
% B = bidiag(A)
% [U,B,V] = bidiag(A)
%
% Computes the bidiagonalization of the m-times-n matrix A with m >= n:
%     A = U*B*V' ,
% where B is an upper bidiagonal n-times-n matrix, and U and V have
% orthogonal columns.  The matrix B is stored as a sparse matrix.

% Reference: L. Elden, "Algorithms for regularization of ill-
% conditioned least-squares problems", BIT 17 (1977), 134-145.

% Per Christian Hansen, IMM, Sept. 13, 2001.

% Initialization.
[m,n] = size(A);
if (m < n), error('Illegal dimensions of A'), end
B = sparse(n,n);
if (nargout> 1), U = [eye(n);zeros(m-n,n)]; betaU = zeros(n,1); end
if (nargout==3), V = eye(n); betaV = zeros(n,1); end

% Bidiagonalization; save Householder quantities.
if (m > n), k_last = n; else k_last = n-1; end
for k=1:k_last

  [B(k,k),beta,A(k:m,k)] = gen_hh(A(k:m,k));
  if (k < n), A(k:m,k+1:n) = app_hh(A(k:m,k+1:n),beta,A(k:m,k)); end
  if (nargout>1), betaU(k) = beta; end

  if (k < n-1)
    [B(k,k+1),beta,v] = gen_hh(A(k,k+1:n).'); A(k,k+1:n) = v.';
    A(k+1:m,k+1:n) = app_hh(A(k+1:m,k+1:n)',beta,A(k,k+1:n)')';
    if (nargout==3), betaV(k) = beta; end
  elseif (k == n-1)
    B(n-1,n) = A(n-1,n);
  end

end

% Save bottom element if A is square.
if (k_last < n), B(n,n) = A(n,n); end

% Compute U if wanted.
if (nargout>1)
  for k=k_last:-1:1
    U(k:m,k:n) = app_hh(U(k:m,k:n),betaU(k),A(k:m,k));
  end
end

% Compute V if wanted.
if (nargout==3)
  for k=n-2:-1:1
    V(k+1:n,k:n) = app_hh(V(k+1:n,k:n),betaV(k),A(k,k+1:n)');
  end
end

if (nargout < 2), U = B; end