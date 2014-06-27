% Inference example for sparse regression example using glm-ie.
% Exercises opts.outerExact.
%
% We use the Nordborg Flowering Time dataset from http://walnut.usc.edu/.
%   X is a binary matrix of 166 individuals x 5k SNPs selected by prior
%     knowledge
%   y is a real vector with flowering times between 21 and 200
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

load regress, y = (y-mean(y))/std(y); [m,n] = size(X);                % get data

s2 = 0.0078;                                                          % variance
B = 1;  tau = 125;  t = 0;                                  % sparsity transform

% Gaussian or non-Gaussian i.e. sparse weight prior
% select by uncommenting the desired one
% pot = @(s) potT(s,2);            % Student's t
% pot = @potLaplace;               % Laplace
pot = @(s) potExpPow(s,1.3);     % Exponential Power family

fprintf('Do sparse regression using a (sparse) '), str = func2str(pot);
fprintf('%s weight prior.\n',str(8:end))

% split into training and test set
itr = randperm(m); ite = itr(1:end/2); itr = itr(end/2+1:end);

%% Inference
opts.innerOutput = 1;
opts.outerOutput = 1;
opts.outerNiter = 3;           % Number of outer loop iterations
opts.outerMethod = 'lanczos';  % Lanczos marginal variance approximation
% opts.outerMethod = 'sample';   % Monte Carlo marginal variance approximation
% opts.outerMethod = 'full';     % Exact marginal variance computation (slow)
% opts.outerMethod = 'woodbury'; % Exact marginal variance computation (fast)

opts.innerMVM =  40;           % number CG steps
opts.innerVBpls = 'plsCG';     % PLS algorithm, also LBFGS if compiled
[m,ga,b,z,zu,nlZ,Q,T] = dli(X(itr,:),y(itr),s2,B,t,pot,tau,opts);

% posterior quantities:
% m = E(u) posterior mean estimate
% z  = var(B*u) marginal variance estimate
% nlZ sequence of negative log marginal likelihoods, nlZ(end) is the last one
% zu = var(u) marginal variance estimate
% Q,T yield cov(u) = V = inv(A) \approx Q'*inv(T)*Q the posterior covariance

%% Estimation
opt.nMVM = 50; opt.output = 1; opt.exactNewt = 1;
uhat = feval(opts.innerVBpls,zeros(n,1),X(itr,:),y(itr),B,t,opt,s2,'penVB',pot,tau);

fprintf('mean squared error (MSE)\n  ')
fprintf('inference  %1.2f, ', sum((X(ite,:)*m-y(ite)).^2)/numel(ite) )
fprintf('estimation %1.2f, ', sum((X(ite,:)*uhat-y(ite)).^2)/numel(ite) )
fprintf('\b\b\n')
