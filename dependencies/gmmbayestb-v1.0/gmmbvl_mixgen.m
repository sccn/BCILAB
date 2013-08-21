function [X,T,L1,L2] = gmmbvl_mixgen(n,m,k,d,c,e)
%gmmbvl_mixgen - Gaussian mixture generator 
%
%[X,T] = gmmbvl_mixgen(n,m,k,d,c,e) 
%  n - size of training set 
%  m - size of test set 
%  k - number of components
%  d - dimension
%  c - separation degree
%  e - maximum eccentricity
%returns
%  X - training set (n x d)
%  T - test set (m x d) 

% Nikos Vlassis, 2000
% for definitions see (Dasgupta, 1999)

%
% $Name:  $

R=zeros(k,d^2);

% mixing weights
while 1
  W = rand(k,1); 
  W = W / sum(W);
  if all(W > 1/(4*k)) 
    break;
  end
end

% create c-separated Gaussian clusters of maximum eccentricity e
trials = 1;
while 1
  X = [];
  T = [];
  M = randn(k,d)*sqrt(k)*sqrt(c)*trials/10;
  Trace = zeros(k,1);
  for j = 1:k
    U = rand(d,d)-0.5; 
    U = sqrtm(inv(U*U')) * U;
    L = diag(rand(d,1)*(e-1)+1).^2/100;
    msg = 1;
    while msg
      [C,msg] = chol(U*L*U');
    end
    R(j,:)=C(:)';

    nj = ceil(n*W(j));
    Xj = randn(nj,d) * C;
    X = [X; repmat(M(j,:),nj,1) + Xj];
    Trace(j) = trace(cov(Xj));

    mj = ceil(m*W(j));
    Tj = randn(mj,d) * C;
    T = [T; repmat(M(j,:),mj,1) + Tj];

  end
  % check degree of separation
  error = 0;
  for i = 1:k-1
    for j = i+1:k
      if norm(M(i,:)-M(j,:)) < c * sqrt(max(Trace(i),Trace(j)))
        error = 1;
      end
    end
  end
  if ~error
    break;
  end
  trials = trials + 1;
end

L = gmmbvl_em_gauss(X,M,R);
F = L*W; 
F(find(F < eps)) = eps;
L1 = mean(log(F)); 
if ~isempty(T) 
  L = gmmbvl_em_gauss(T,M,R);
  F = L*W; 
  F(find(F < eps)) = eps;
  L2 = mean(log(F));
end
