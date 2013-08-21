function [W,M,R,P,sigma] = gmmbvl_em_init_km(X,k,dyn)
%gmmbvl_em_init_km - initialization of EM for Gaussian mixtures 
%
%[W,M,R,P,sigma] = gmmbvl_em_init_km(X,k,dyn)
%  X - (n x d) matrix of input data 
%  k - initial number of Gaussian components
%  dyn - if 1 then perform dynamic component allocation else normal EM 
%returns
%  W - (k x 1) vector of mixing weights
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of Cholesky submatrices of components covariances
%  P - (n x k) the posteriors to be used in EM step after initialization
%  of priors, means, and components covariance matrices

% Nikos Vlassis & Sjaak Verbeek 2002

%
% $Name:  $

[n,d] = size(X);

[tmp,M,tmp2] = gmmbvl_kmeans(X,[],k,0,0,0,0);
[D,I]        = min(gmmbvl_sqdist(M',X'),[],1);

% mixing weights
W = zeros(k,1);
for i=1:k
	W(i) = length(find(I==i))/n;
end

% covariance matrices
R = zeros(k,d^2);
if k > 1
	for j = 1:k
		J = find(I==j);
		if length(J)>2*d;
			Sj = cov(X(J,:));
		else
			Sj = cov(X);
		end
		Rj = chol(Sj);
		R(j,:) = Rj(:)';
	end
else
	S = cov(X);
	R = chol(S);
	R = R(:)';
end

% compute likelihoods L (n x k)
L = gmmbvl_em_gauss(X,M,R);

% compute mixture likelihoods F (n x 1)
F = L * W;
F(find(F < eps)) = eps;

% compute posteriors P (n x k)
P = L .* repmat(W',n,1)  ./ repmat(F,1,k);

sigma = 0.5 * (4/(d+2)/n)^(1/(d+4)) * sqrt(norm(cov(X)));
