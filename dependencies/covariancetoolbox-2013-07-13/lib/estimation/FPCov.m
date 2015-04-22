function C = FPCov(X)
%Fixed Point covariance matrix estimation
N = size(X,1);
M = size(X,2);
C = NormalizedSCM(X);
critere = 1 ;

Max_iter = 50;
iter = 0;
while ((critere)>10^-3)&&(iter<Max_iter)
    iter = iter +1;
    C_precedant = C;
    X2 = X./sqrt(repmat(diag(X*inv(C_precedant)*X'),1,M));
    C = (M/N)*(X2')*X2;
    critere = norm(C-C_precedant,'fro');
end

