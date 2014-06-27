% Inference example for logistic regression using glm-ie demonstrating the 
% sparsity occuring in estimation but similar classification performance in
% inference where probabilities can readily be evaluated.
%
% We use the well known a9a dataset from the UCI machine learning repository
% with 32561 instances of dimension 123.
% B data matrix
% c binary labels
%
% (c) by Hannes Nickisch, Philips Research, 2013 August 30
clear all, close all

potS = @potLaplace; % Laplace weight prior

% load and split into train and test set
load classify                               % load a part of the UCI a9a dataset
nte = 30000;                        % say how many test examples we wich to keep
Bte = B(1:nte,:); B = B(nte+1:end,:); cte = c(1:nte); c = c(nte+1:end);  % split
[q,n] = size(B);

% nearest neighbour classification
dp = sum(bsxfun(@minus, Bte, full(mean(B(c==+1,:)))).^2,2);  % dist. to pos. ctr
dm = sum(bsxfun(@minus, Bte, full(mean(B(c==-1,:)))).^2,2);  % dist. to neg. ctr
cc = 2*double(dp<dm)-1;
fprintf('nearest neighbor accuracy=%1.2f%%\n',100*sum(cte==cc)/numel(cte))

% classification
fprintf('Do regression using a ')
str = func2str(potS);
fprintf('%s',str(8:end))
m = 0; X = zeros(0,n); y = zeros(0,1); s2 = 1;         % no Gaussian part at all
tau  = [2*ones(n,1); 5*c.*ones(q,1)];
B = [eye(n); B]; t = 0;
pot = @(s,varargin) potCat(s, varargin{:}, {potS,@potLogistic}, {1:n,n+(1:q)} );

fprintf(' weight prior.\n')

%% Inference
opts.innerOutput = 1;
opts.outerOutput = 1;
opts.outerNiter  = 3;         % Number of outer loop iterations
opts.outerMVM   = 123;        % Lanczos steps
opts.outerMethod = 'lanczos';
opts.innerMVM   =  50;        % CG steps
opts.innerVBpls = 'plsTN';    % PLS algorithm, use LBFGS if compiled
opts.innerType = 'VB';
% opts.innerType = 'EP'; opts.innerEPeta = 0.9;
[uinf,ga,b,z,zu,nlZ] = dli(X,y,s2,B,t,pot,tau,opts);
sdinf = sqrt(zu);

%% Estimation
opt.nMVM = 150; opt.output = 1;
u0 = zeros(n,1);
uest = feval(opts.innerVBpls,u0,X,y,B,t,opt,s2,'penVB',pot,tau/sqrt(s2));

fprintf('accuracies on test set\n')
acc_inference  = 100*sum(cte==sign(Bte*uinf))/numel(cte)
acc_estimation = 100*sum(cte==sign(Bte*uest))/numel(cte)

% ci-f1.png
sz = [500,300]; figure('Position',[50,50,sz]), set(gca,'FontSize',16)
plot([-1,1],[-1,1],'r'), hold on
for i=1:length(uinf)
  plot(uinf(i)+sdinf*[-1,1],uest(i)*[1,1],'c')
end
plot(uinf,uest,'bo')
plot(uinf,uest,'k.')
axis([-1,1,-1,1])
xlabel('u_{VB/EP}'), ylabel('u_{MAP}')
title('sparsity pattern: inference vs. estimation')