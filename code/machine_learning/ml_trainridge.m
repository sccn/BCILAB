function model = ml_trainridge(varargin)
% Perform ridge regression.
% Model = ml_trainridge(Trials, Targets, Lambda, Options...)
%
% Ridge regression [1] is a generalization of basic multivariate regression which in addition
% includes a regularization parameter (lambda), which, when chosen as > 0 governs the amount of
% shrinkage applied to the solution. This is an effective way to carefully contrain the model
% complexity in order to not overfit spurious noise in the training data. It can also be interpreted
% as a Gaussian prior on the model weights when ridge regression is viewed as a Bayesian Maximum a
% Posteriori (MAP) estimator. This type of regularization is also called Tikhonov regularization.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%                  may also contain weights, in the form {Targets,Weights}; supported by l1 & l2 variants
%
%   Lambda       : regularization parameter (default: 1)
%                  a comprehensive search interval would be 2.^(-5:2:15)
%
%   Options  : optional name-value parameters to control the training details:
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std')
%
% Out:
%   Model   : a linear model;
%             w is the linear weights, including the bias as last element;
%
% See also:
%   ml_predictridge
%
% Examples:
%   % learn a ridge regression model using defaults
%   model = ml_trainlogreg(trials,targets)
%
%   % as before, but using a non-default regularization parameter (here: 0.1)
%   model = ml_trainlogreg(trials,targets,0.1)
%
%   % as before, but also use a non-default scaling option (center the data)
%   model = ml_trainlogreg(trials,targets,0.1,'scaling','center')
%
%   % use ml_trainridge with the parameter search function, and search over a range of possible lambda values
%   model = utl_searchmodel({trials,targets},'args',{{'ridge',search(2.^(-5:1:10)),'scaling','center'}})
%
% References:
%   [1] Hastie, T., Tibshirani, R., & Friedman, J. 
%       "The Elements of Statistical Learning: Data Mining, Inference, and Prediction" 
%       Springer Series in Statistics (2009)
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2013-11-17

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'lambda','Lambda'}, 1, [0 2^-7 2^15 Inf], 'Regularization parameter. Reasonable range: 2.^(-5:2:15), greater is stronger.'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.','cat','Miscellaneous'));

% scale the data
sc_info = hlp_findscaling(trials,scaling); %#ok<*NODEF>
trials = hlp_applyscaling(trials,sc_info);

trials = [trials ones(size(trials,1),1)];
model.w = (trials'*trials + lambda*eye(size(trials,2)))\(trials'*targets); 

model.sc_info = sc_info;
