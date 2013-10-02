function X = robust_mean(Y,rho,iters)
% Perform a robust mean under the Huber loss function.
% x = robust_mean(Y,rho,iters)
%
% Input:
%   Y : MxN matrix over which to average (columnwise)
%   rho : augmented Lagrangian variable (default: 1)
%   iters : number of iterations to perform (default: 1000)
%
% Output:
%   x : 1xN vector that is the roust mean of Y
%
% Based on the ADMM Matlab codes also found at:
%   http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-09-26

if ~exist('rho','var')
    rho = 1; end
if ~exist('iters','var')
    iters = 1000; end

m = size(Y,1);
if m==1
    X = Y;
else
    mu = sum(Y)/m;
    Z = zeros(size(Y)); U = Z;
    for k = 1:iters
        X = mu + sum(Z - U)/m;
        D = bsxfun(@minus,X,Y - U);
        Z = (rho/(1+rho) + (1/(1+rho))*max(0,(1-(1+1/rho)./abs(D)))).*D;
        U = D - Z;
    end
end