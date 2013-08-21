function model = ml_traingauss(varargin)
% Learn a predictive model via a robust Gaussian Bayes classifier (with feature selection).
% Model = ml_trainlda(Trials, Targets, Lambda, Options...)
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train; 
%                  can include a per-trial weighting (if a cell array {Targets,Weights})
%
%   Options...   : optional name-value pairs, with possible names:
%                  'dimensions': if given, reduce data dimensionality to this number of dimensions using a strong feature selector (default: [])
%                  'blend' : covariance blending parameter, blends covariance matrices (default: 1 = all covs identical)
%                            can be selected by a parameter search, e.g., by specifying search(0:0.1:1)
%                  'normprobs': produce normalized probabilities (0/1) (default: 1)
%                  'scaling': data pre-scaling (default: 'center')
%
% Out:
%   Model   : a predictive mdoel
%
% Examples:
%   % learn a simple Gaussian Bayes classifier
%   model = ml_traingauss(trials,targets)
%
%   % as before, but reduce the dimensions of the data to 10
%   model = ml_traingauss(trials,targets,'dimensions',10)
%
%   % as before, but allow for non-identical covariance matrices
%   model = ml_traingauss(trials,targets,'dimensions',10,'blend',0.5)
%
%   % like before, but search over the dimensions and covariance blending parameter using parameter search
%   model = utl_searchmodel({trials,targets},'args',{{'gauss','dimensions',search(2:2:20),'blend',search(0:0.1:1)}})
%   
% See also:
%   ml_predictgauss
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03

arg_define([0 2],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'dimensions','Dimensions'}, [], [], 'Reduce data dimentions to n. If given, reduce data dimensionality to this number of dimensions using a strong feature selector.','shape','scalar'), ...
    arg({'blend','CovBlending'}, 1, [], 'Degree of blending between per-class covariances. Blends covariance matrices (0=all covs independent, 1=all covs identical) and can be selected by a parameter search, e.g., by specifying search(0:0.1:1)'), ...
    arg({'normprobs','NormalizedProbs'}, true, [], 'Produce normalized probabilities. If turned off, the classifier can be used for density estimation (e.g. novelty discovery).'), ...
    arg({'scaling','Scaling'}, 'center', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg({'startcost','InitialCost'},2,[],'Initial cost for dimensionality reduction. If dimensionality reduction is requested, this is the initial cost which will be lowered until the desired number of dimensions is reached.'));

% scale the data
model.sc_info = hlp_findscaling(trials,scaling);
trials = hlp_applyscaling(trials,model.sc_info);

% select features
if ~isempty(dimensions)
    % start with some fairly high cost
    cost = startcost;
    while 1
        model.selector = ml_trainlogreg(trials,targets,cost,'variant','l1','scaling','');
        w = sum(abs(model.selector.w(:,1:end-1)));
        % reduce cost until we have at last #dimensions non-zero weights
        if nnz(w) >= dimensions
            break; end
        cost = cost/2;        
    end
    [dummy,I] = sort(w,'descend'); %#ok<ASGLU>
    model.feature_mask = full(sparse(1,I(1:dimensions),true,1,length(I)));
    model.feature_weight = w;
    trials = trials(:,model.feature_mask);
else
    model.feature_mask = true(1,size(trials,2));
    model.feature_weight = ones(1,size(trials,2));
end

% obtain weights
if iscell(targets)
    [targets,weights] = deal(targets{:});
else
    weights = ones(size(targets,1),1);
end

% estimate the distribution of each class, using a weighted estimator (with shrinkage)...
model.classes = unique(targets);
if length(model.classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.'); end
for c=1:length(model.classes)
    w = weights(targets == model.classes(c));
    X = trials(targets == model.classes(c),:);
    model.mu{c} = mean_w(X,w);
    model.sigma{c} = cov_shrink(X,w);
end

% blend covariance matrices
if blend ~= 0
    meansigma = mean(cat(3,model.sigma{:}),3);
    for c=1:length(model.classes)
        model.sigma{c} = (1-blend)*model.sigma{c} + blend*meansigma; end
end

model.normprobs = normprobs;
