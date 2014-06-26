function Kt = center_gram_matrix(K);
n = size(K,1);
sumK = sum(K,2);
Kt = K - repmat(sumK,1,n)/n - repmat(sumK,1,n)'/n + sum(sumK)/n/n;
Kt = .5 * ( Kt + Kt');