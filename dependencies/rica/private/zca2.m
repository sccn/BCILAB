function [Z,V,E,D] = zca2(X)
eps = 1e-4;
% Calculate the eigenvalues and eigenvectors of the new covariance matrix.
covarianceMatrix = X*X'/size(X,2);
[E, D] = eig(covarianceMatrix);

% Sort the eigenvalues  and recompute matrices
[~,order] = sort(diag(-D));
E = E(:,order);
d = diag(D); 
dsqrtinv = real((d + ones(size(d))*eps).^(-0.5));
Dsqrtinv = diag(dsqrtinv(order));
D = diag(d(order));
V = E*Dsqrtinv*E';
Z = V*X;

