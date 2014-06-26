function [W,S,m] = runSoftICA(patches,numcomps,chanlocs)


%% Whiten and rescale the data
sphere = inv(real(sqrtm(cov_robust(patches'))));
patches = sphere*patches;
m = sqrt(sum(patches.^2) + (1e-8));
x = bsxfunwrap(@rdivide,patches,m);

%% Run the optimization

params.lambda = 0.05;
params.numFeatures = numcomps;
params.epsilon = 1e-5;
params.n = size(patches,1);

%configure minFunc
options.Method = 'lbfgs';
options.MaxFunEvals = Inf;
options.MaxIter = 300;
%options.display = 'off';
%options.outputFcn = 'showBases';

% initialize with random weights
randTheta = randn(params.numFeatures,params.n)*0.01;  % 1/sqrt(params.n);
randTheta = randTheta ./ repmat(sqrt(sum(randTheta.^2,2)), 1, size(randTheta,2)); 
randTheta = randTheta(:);

% optimize
[opttheta, cost, exitflag] = minFunc( @(theta) softICACost(theta, x, params), randTheta, options);   % Use x or xw 

% display result
W = reshape(opttheta,params.numFeatures,params.n);
S = sphere;

