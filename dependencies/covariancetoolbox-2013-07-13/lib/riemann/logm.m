function Out = logm(X)

[V D] = eig(X);
Out = V*diag(log(diag(D)))*V';
