function Out   = RiemannLogMap(P,X)

[U Delta] = eig(P);
G = U*sqrt(Delta);
Y = inv(G)*X*inv(G)';
[V Sigma] = eig(Y);
Out = (G*V)*diag(log(diag(Sigma)))*(G*V)';
