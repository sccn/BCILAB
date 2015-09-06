function [x_k,rho,eta] = ttls(V1,k,s1)
%TTLS Truncated TLS regularization.
%
% [x_k,rho,eta] = ttls(V1,k,s1)
%
% Computes the truncated TLS solution
%    x_k = - V1(1:n,k+1:n+1)*pinv(V1(n+1,k+1:n+1))
% where V1 is the right singular matrix in the SVD of the matrix
%    [A,b] = U1*diag(s1)*V1' .
%
% If k is a vector, then x_k is a matrix such that
%    x_k = [ x_k(1), x_k(2), ... ] .
% If k is not specified, k = n is used.
%
% The solution norms and TLS residual norms corresponding to x_k are
% returned in eta and rho, respectively.  Notice that the singular
% values s1 are required to compute rho.

% Reference: R. D. Fierro, G. H. Golub, P. C. Hansen and D. P. O'Leary,
% "Regularization by truncated total least squares", SIAM J. Sci. Comput.
% 18 (1997), 1223-1241.

% Per Christian Hansen, IMM, 03/18/93.

% Initialization.
[n1,m1] = size(V1); n = n1-1;
if (m1 ~= n1), error('The matrix V1 must be square'), end
if (nargin == 1), k = n; end
lk = length(k);
if (min(k) < 1 | max(k) > n)
  error('Illegal truncation parameter k')
end
x_k = zeros(n,lk);
if (nargout > 1)
  if (nargin < 3)
    error('The singular values must also be specified')
  end
  ns = length(s1); rho = zeros(lk,1);
end
if (nargout==3), eta = zeros(lk,1); end

% Treat each k separately.
for j=1:lk
  i = k(j);
  v = V1(n1,i+1:n1); gamma = 1/(v*v');
  x_k(:,j) = - V1(1:n,i+1:n1)*v'*gamma;
  if (nargout > 1), rho(j) = norm(s1(i+1:ns)); end
  if (nargout == 3), eta(j) = sqrt(gamma - 1); end
end