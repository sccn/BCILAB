function C = NormalizedSCM(X)
N = size(X,1);
M = size(X,2);
X2 = X./sqrt(repmat(diag(X*X'),1,M));
C = (M/N)*(X2')*X2;