function [V,E] = quadinv(N,NW)

% calculates the quadratic inverse eigenvectors
%
% Input
% N: number of samples
% NW: time-bandwidth product
%
% Output
% V: The quadratic inverse eigenvectors
% E: the quadratic inverse eigenvalues


x = 0:N-1;
indx = find(x);
y = ones(size(x))*(2*NW/N)^2;

y(indx) = (sin(2*pi*NW*x(indx)/N)./(pi*x(indx))).^2;

M = toeplitz(y);
[C,D] = eig(N*M);

tmp = diag(D);

K = 4*NW;

V = sqrt(N)*C(:,N:-1:N-K+1);

E = tmp(N:-1:N-K+1);

for ii =1:K
  if(sum((N-1-2*x)'.*V(:,ii))<0)
    V(:,ii) = flipud(V(:,ii));
  end;
end;

return;
