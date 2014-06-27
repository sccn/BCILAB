% Estimation example for logistic regression using glm-ie.
%
% We use the well known a9a dataset from the UCI machine learning repository
% with 32561 instances of dimension 123.
% B data matrix
% c binary labels
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

% Gaussian or non-Gaussian i.e. sparse weight prior
% select by uncommenting the desired one
% gauss = 1;                                   % Gaussian
gauss = 0; potS = @(s) potT(s,2);            % Student's t
% gauss = 0; potS = @(s) potLaplace(s);        % Laplace
% gauss = 0; potS = @(s) potExpPower(s,0.7);   % Exponential Power family

% load and split into train and test set
load classify                               % load a part of the UCI a9a dataset
nte = 20000;                        % say how many test examples we wich to keep
Bte = B(1:nte,:); B = B(nte+1:end,:); cte = c(1:nte); c = c(nte+1:end);  % split
[q,n] = size(B); t = 0;

% nearest neighbour classification
dp = sum(bsxfun(@minus, Bte, full(mean(B(c==+1,:)))).^2,2);  % dist. to pos. ctr
dm = sum(bsxfun(@minus, Bte, full(mean(B(c==-1,:)))).^2,2);  % dist. to neg. ctr
cc = 2*double(dp<dm)-1;
fprintf('nearest neighbor accuracy=%1.2f%%\n',100*sum(cte==cc)/numel(cte))

% classification
fprintf('Do classification using a ')
if gauss
  fprintf('Gaussian')
  X = eye(n); y = zeros(n,1); m = n; s2 = 1;
  tau = c; pot = @(s) potLogistic(s);
else
  str = func2str(potS);
  fprintf('%s',str(8:end))
  m = 0; X = zeros(0,n); y = zeros(0,1); s2 = 1;       % no Gaussian part at all
  tau  = [2*ones(n,1); 5*c.*ones(q,1)];
  B = [eye(n); B];
  pot = @(s) [potS(s(1:n)); potLogistic(s(n+(1:q)))];
end
fprintf(' weight prior.\n')

% optimisation parameters
opt.nMVM = 250; opt.output = 1;
u0 = zeros(n,1);

%% Estimation: apply 3 different PLS schemes
plsList = {'TN','CGBT','CG'}; % also LBFGS
pars = {u0,X,y,B,t,opt,s2,'penVB',pot,tau};
for i=1:length(plsList)
  fprintf('PLS optimisation using %s.\n',plsList{i})
  tic, [u{i},phi(i)] = feval(['pls',plsList{i}],pars{:}); tt(i) = toc;
  cc = sign(Bte*u{i}); acc(i) = 100*sum(cte==cc)/numel(cte);
end

fprintf('accuracies in %%\n  ')
for i=1:length(plsList), fprintf('%s %1.2f%%, ',plsList{i},acc(i)), end
fprintf('\b\b\n')

fprintf('running times\n  ')
for i=1:length(plsList), fprintf('%s %1.2fs, ',plsList{i},tt(i)), end
fprintf('\b\b\n')

fprintf('objective function values\n  ')
for i=1:length(plsList), fprintf('%s %1.4e, ',plsList{i},phi(i)), end
fprintf('\b\b\n')