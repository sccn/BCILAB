% demo script for regression with HKL
clear all


% fixing the seed of the random generators
seed=1;
randn('state',seed);
rand('state',seed);

% toy example characteristics
p = 1024;           % total number of variables (used to generate a Wishart distribution)
psub = 8;          % kept number of variables = dimension of the problem
n = 512;            % number of observations
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

disp('Learning with the Hermite kernels');

% HKL with hermite polynomial decomposition
[outputs,model,accuracies] = hkl(X(trainset,:),Y(trainset),lambdas,'square','hermite',[ .5 3 .1 4 ],...
	'Xtest',X(testset,:),'Ytest',Y(testset),'maxactive',300,'memory_cache',1e9,'display',0);

% perform testing on the test set provided separately (checking that predictions are the same, and showing how to do so)
accuracies_test = hkl_test(model,outputs,X(testset,:),'Ytest',Y(testset));


subplot(1,2,1);


plot(-log10(lambdas),accuracies.testing_error,'r'); hold on;
plot(-log10(lambdas),accuracies.training_error,'b'); hold on;
legend('test','train');
xlabel('log_{10}(\lambda)');
title('hermite kernels');


disp('Learning with the ANOVA kernels');

% HKL with ANOVA kernel
[outputs,model,accuracies] = hkl(X(trainset,:),Y(trainset),lambdas,'square','anova',[ .0625 .1 8 30],...
	'Xtest',X(testset,:),'Ytest',Y(testset),'maxactive',300,'memory_cache',1e9);


disp('Learning with the ANOVA kernels (given as kernel matrices)');
% does the same where the kernel matrices are given directly (for checking that the models are the same and showing how to do so)
[n , p ] = size(X);
temppmean = mean(X(trainset,:),1);	% stored to allow testing on new data points
Xtemp = X - repmat(temppmean,n,1);
tempstd = sqrt( mean(Xtemp(trainset,:).^2,1) );
Xtemp = Xtemp ./ repmat(tempstd,n,1);
Ks = zeros(ntrain*(ntrain+1)/2,size(X,2),1,'single');
Kstest = zeros(ntrain*ntest,size(X,2),1,'single');
b = .0625;
for i=1:size(X,2)
	Ks(:,i,1) = symmetric_vectorize_single( single( exp( - b * sqdist(Xtemp(trainset,i)',Xtemp(trainset,i)' ) ) ));
	temp =  single(exp( - b * sqdist(Xtemp(testset,i)',Xtemp(trainset,i)' ) ) );
	Kstest(:,i,1) = temp(:);
end



[outputs_direct,model_direct,accuracies_direct] = hkl(Ks,Y(trainset),lambdas,'square','base kernels',[ size(X,2) 1 .1 8 30],...
	'Xtest',Kstest,'Ytest',Y(testset),'maxactive',300,'memory_cache',1e9);

accuracies_direct_test = hkl_test(model_direct,outputs_direct,Kstest,'Ytest',Y(testset));


subplot(1,2,2);
plot(-log10(lambdas),accuracies.testing_error,'r'); hold on;
plot(-log10(lambdas),accuracies.training_error,'b'); hold on;
legend('test','train');
xlabel('log_{10}(\lambda)');
title('anova kernels');






