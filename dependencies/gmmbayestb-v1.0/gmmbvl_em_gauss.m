function L = gmmbvl_em_gauss(X,M,R)
%gmmbvl_em_gauss - compute likelihoods for all points and all components
%
%L = gmmbvl_em_gauss(X,M,R)
%  X - (n x d) matrix of input data
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of Cholesky submatrices of components covariances
%      in vector reshaped format. To get the covariance of component k:
%      Rk = reshape(R(k,:),d,d); S = Rk'*Rk;
%returns 
%  L - (n x k) likelihoods of points x_n belonging to component k

% Nikos Vlassis, 2000

%
% $Name:  $

[n,d] = size(X);
k = size(M,1);

L = zeros(n,k); 
for j = 1:k 

  % Cholesky triangular matrix of component's covariance matrix
  Rj = reshape(R(j,:),d,d);        
  
  % We need to compute the Mahalanobis distances between all inputs
  % and the mean of component j; using the Cholesky form of covariances
  % this becomes the Euclidean norm of some new vectors 
  New = (X - repmat(M(j,:),n,1)) * inv(Rj);
  Mah = sum(New.^2,2);

  L(:,j) = (2*pi)^(-d/2) / det(Rj) * exp(-0.5*Mah);
end
