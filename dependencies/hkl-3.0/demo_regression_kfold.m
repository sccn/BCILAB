% demo script for regression with HKL
clear all


% fixing the seed of the random generators
seed=1;
randn('state',seed);
rand('state',seed);

% toy example characteristics
p = 1024;           % total number of variables (used to generate a Wishart distribution)
psub = 32;          % kept number of variables = dimension of the problem
n = 256;            % number of observations
s = 4;              % number of relevant variables
noise_std = .2;		% standard deviation of noise
proptrain = .5;     % proportion of data kept for training (the rest is used for testing)


% generate random covariance matrix from a Wishart distribution
Sigma_sqrt = randn(p,p);
Sigma = Sigma_sqrt' * Sigma_sqrt;


% normalize to unit trace and sample
diagonal = diag(Sigma);
Sigma = diag( 1./diagonal.^.5) * Sigma * diag( 1./diagonal.^.5);
Sigma_sqrt =   Sigma_sqrt * diag( 1./diagonal.^.5);
X = randn(n,p) * Sigma_sqrt;

X = X(:,1:psub);
p=psub;

% generate nonlinear function of X as the sum of all cross-products
J =  1:s;    % select the first s variables
Y = zeros(n,1);
for i=1:s
	for j=1:i-1
		Y = Y + X(:,J(i)) .* X(:,J(j));
	end
end
% normalize to unit standard deviation
Y = Y / std(Y);

% add some noise with known standard deviation
Y =  Y + randn(n,1) * noise_std;


% split data in two groups
ntrain = round(n*proptrain);
ntest = n - ntrain;
rp = randperm(n);
trainset = rp(1:ntrain);
testset  = rp(ntrain+1:end);

% select some regularization parameterst lambdas (from large to small)
lambdas = 10.^[2:-.5:-8];


% HKL with hermite polynomial decomposition
[outputs,model,accuracies,accuracies_kfold,bestlambda_same] = hkl_kfold(10,'same',X(trainset,:),Y(trainset),lambdas,'square','hermite',[ .5 3 .1 4 ],...
	'Xtest',X(testset,:),'Ytest',Y(testset),'maxactive',300,'memory_cache',2e9,'display',0);


estimated_best_error_same = accuracies.testing_error(bestlambda_same)
best_error_same = min(accuracies.testing_error)



[outputs,model,accuracies,accuracies_kfold,bestlambda_scaled] = hkl_kfold(10,'scaled',X(trainset,:),Y(trainset),lambdas,'square','hermite',[ .5 3 .1 4 ],...
	'Xtest',X(testset,:),'Ytest',Y(testset),'maxactive',300,'memory_cache',2e9,'display',0);


estimated_best_error_scaled = accuracies.testing_error(bestlambda_scaled)
estimated_best_error_overfit = min(accuracies.testing_error)



