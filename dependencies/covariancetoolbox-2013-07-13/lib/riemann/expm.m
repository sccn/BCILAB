function Out = expm(X)

[V D] = eig(X);
Out = V*diag(exp(diag(D)))*V';
