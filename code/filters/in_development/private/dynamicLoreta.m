function [J,sigma2,tau2,T,history] = dynamicLoreta(Y, Ut, s2,iLV,sigma2,tau2, options)

%[J,varargout] = dynamicLoreta(V,varargin)
%
% Computes the posterior distribution of the parameters J given some data V. 
% The program solves levels of inference: 1) optimization of parameters J, and
% 2) optimization of hyperparameters sigma and tau. See references for details.
%
% Ut,s2, and iLV are defined as follows: 
%     Y: Nsensors x time points data matrix
%     K: N x P predictor matrix
%     L: sparse P x P square root of the precision matrix 
%     [U,s,V] = svd( K*\L )
%     iLV = inv(L)*V
%     s2  = s.^2
%
% sigma, tau: hyperparameters
% J: current source density (estimated parameters)
% 
%                     p(V|J,sigma)*P(J|tau)
% p(J|V,sigma,tau) = ---------------------- 
%                        p(V|sigma,tau)
% 
%                     /      
% l(sigma, tau) = log | p(V|J,sigma) *p(J|tau)
%                     /
% References:
%   Trujillo-Barreto, N., Aubert-Vazquez, E., Valdes-Sosa, P.A., 2004.
%       Bayesian model averaging in EEG/MEG imaging. NeuroImage 21, 1300???1319
%
%   Yamashita, O., Galka,A., Ozaki, T., Biscay, R. Valdes-Sosa, P.A., 2004.
%       Recursive Penalized Least Squares Solution for Dynamical Inverse Problems
%       of EEG Generation. Human Brain Mapping 21:221–235
%
% Author: Alejandro Ojeda, Syntrogi Inc., Jan-2014

if nargin < 4, error('Not enough input arguments.');end
if nargin < 5, sigma2 = [];end
if nargin < 6, tau2 = [];end
if nargin < 7,
    options = struct('maxTol',1e-3,'maxIter',100,'gridSize',100,'verbose',true,'history',true,'useGPU',false,'initNoiseFactor',0.001);
end
[history.sigma2, history.tau2, history.gcv, history.error] = deal(nan(1,options.maxIter));
error_win = 3;
Y2 = Y.^2;
s = s2.^(0.5);
n = numel(Y);
p = numel(Y);

% Initialize hyperparameters
if isempty(sigma2) || isempty(tau2)
    UtY = Ut*Y;
    tol = max([n p])*eps(max(s));
    lambda2 = logspace(log10(tol),log10(max(s)),options.gridSize);
    gcv = zeros(options.gridSize,1);
    for k=1:options.gridSize
        d = lambda2(k)./(s2+lambda2(k));
        f = mean(diag(d)*UtY,2);
        gcv(k) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), loc = 1;end
    lambda2 = lambda2(loc(end));
    sigma2  = options.initNoiseFactor*(Y(:)'*Y(:))/n;
    tau2    = sigma2*lambda2;
    gcv     = gcv(loc(end));
else
    lambda2 = tau2/sigma2;
    gcv = Inf;
end
history.sigma2(1) = sigma2;
history.tau2(1) = tau2;
history.gcv(1) = gcv;
history.error(1) = inf;


if options.verbose
    fprintf('Iter\tSigma2\t\tLambda2\t\tDf\t\tHyperp. Error\tGCV\n');
end
for it=2:options.maxIter
        
    % Computing hat matrix and mse
    H   = Ut'*diag(s2./(s2+lambda2))*Ut;
    mse = mean(sum((Y - H*Y).^2));
        
    % Computing GCV
    gcv = mse/(1-trace(H)/n)^2;
    history.gcv(it) = gcv;   
    
    % Updating hyperparameters
    sigma2  = estimateAlpha(Y2,s2,lambda2);
    lambda2 = estimateLambda2(Y2,sigma2);
    
    history.sigma2(it) = sigma2;
    history.tau2(it)   = sigma2/lambda2;
    if it-error_win < 1
        err = 0.5*std(history.sigma2(1:it)) + 0.5*std(history.tau2(1:it));
    else
        err = 0.5*std(history.sigma2(it-error_win:it)) + 0.5*std(history.tau2(it-error_win:it));
    end
    history.error(it) = err;
    
    if options.verbose
        %disp([num2str(it-1) ' => sigma2: ' num2str(sigma2) '  lambda2: ' num2str(lambda2) ' df: ' num2str( (1-trace(H)/n) ) ' hyperp. error: ' num2str(err) ' gcv: ' num2str(gcv)]);
        fprintf('%i\t%e\t%e\t%e\t%e\t%e\n',it-1,sigma2,lambda2,1-trace(H)/n,err,gcv);
    end
    if err < options.maxTol, break;end
end
if it == options.maxIter, warning('Maximum iteration reached. Failed to converge.');end
if options.verbose
    fprintf('\n')
end
if ~options.history
    history = [];
end
tau2 = sigma2/lambda2;

% parameters's estimation
T = iLV*diag(s./(s2+lambda2))*Ut;
J = T*Y;
end

%--
function sigma2 = estimateAlpha(Y2,s2,lambda2)
sigma2 = mean(mean(bsxfun(@times,Y2,lambda2./(s2+lambda2))));
end
%--
function lambda2 = estimateLambda2(Y2,sigma2)
lambda2 = sigma2/mean(mean(Y2));
end

%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
end