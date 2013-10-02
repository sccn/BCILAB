function x = robust_fit(A,y,rho,iters)
% Perform a robust linear regression using the Huber loss function.
% x = robust_fit(A,y,rho,iters)
%
% Input:
%   A : design matrix
%   y : target variable
%   rho : augmented Lagrangian variable (default: 1)
%   iters : number of iterations to perform (default: 1000)
%
% Output:
%   x : solution for x
%
% Notes:
%   solves the following problem via ADMM for x:
%     minimize 1/2*sum(huber(A*x - y))
%
% Based on the ADMM Matlab codes also found at:
%   http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-03-04

if ~exist('rho','var')
    rho = 1; end
if ~exist('iters','var')
    iters = 1000; end

Aty = A'*y;
L = sparse(chol(A'*A,'lower')); U = L';
z = zeros(size(y)); u = z;
for k = 1:iters
    x = U \ (L \ (Aty + A'*(z - u)));
    d = A*x - y + u;
    z = rho/(1+rho)*d + 1/(1+rho)*max(0,(1-(1+1/rho)./abs(d))).*d;
    u = d - z;
end
