function [L,W] = get_l(n,d)
%GET_L Compute discrete derivative operators.
%
% [L,W] = get_l(n,d)
%
% Computes the discrete approximation L to the derivative operator
% of order d on a regular grid with n points, i.e. L is (n-d)-by-n.
%
% L is stored as a sparse matrix.
%
% Also computes W, an orthonormal basis for the null space of L.

% Per Christian Hansen, IMM, 02/05/98.

% Initialization.
if (d<0), error ('Order d must be nonnegative'), end

% Zero'th derivative.
if (d==0), L = speye(n); W = zeros(n,0); return, end

% Compute L.
c = [-1,1,zeros(1,d-1)];
nd = n-d;
for i=2:d, c = [0,c(1:d)] - [c(1:d),0]; end
L = sparse(nd,n);
for i=1:d+1
  L = L + sparse(1:nd,(1:nd)+i-1,c(i)*ones(1,nd),nd,n);
end

% If required, compute the null vectors W via modified Gram-Schmidt.
if (nargout==2)
  W = zeros(n,d);
  W(:,1) = ones(n,1);
  for i=2:d, W(:,i) = W(:,i-1).*(1:n)'; end
  for k=1:d
     W(:,k) = W(:,k)/norm(W(:,k));
     W(:,k+1:d) = W(:,k+1:d) - W(:,k)*(W(:,k)'*W(:,k+1:d));
  end
end