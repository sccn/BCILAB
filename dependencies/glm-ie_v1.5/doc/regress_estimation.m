% Estimation example for sparse regression example using glm-ie.
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

% penalty
pen = 'penAbsSmooth';                 % equivalent to {'penVB','potLaplace',tau}
lam = s2*tau;
fprintf('Do sparse regression\n')

% optimisation parameters
opt.nMVM = 25;                         % number of matrix vector multiplications
opt.output = 1;                       % whether we show output of the optimisers
opt.exactNewt = 1;         % shall the Newton step in plsTN be computed exactly?
u0 = zeros(n,1);                                                % starting value

% split into training and test set
itr = randperm(m); ite = itr(1:end/2); itr = itr(end/2+1:end);

% apply 4 different PLS estimation schemes and predict
plsList = {'TN','CGBT','CG','BB'}; % also LBFGS
for i=1:length(plsList)
  fprintf('PLS optimisation using %s.\n',plsList{i})
  pls = ['pls',plsList{i}];
  tic, [u{i},phi(i)] = feval(pls,u0,X(itr,:),y(itr),B,t,opt,lam,pen);
  tt(i) = toc;
  mse(i) = sum((X(ite,:)*u{i}-y(ite)).^2)/numel(ite);
end

fprintf('mean squared error (MSE)\n  ')
for i=1:length(plsList), fprintf('%s %1.3e, ',plsList{i},mse(i)), end
fprintf('\b\b\n')

fprintf('running times\n  ')
for i=1:length(plsList), fprintf('%s %1.2fs, ',plsList{i},tt(i)), end
fprintf('\b\b\n')

fprintf('objective function values\n  ')
for i=1:length(plsList), fprintf('%s %1.4e, ',plsList{i},phi(i)), end
fprintf('\b\b\n')