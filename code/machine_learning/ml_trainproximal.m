function model = ml_trainproximal(varargin)
% Learn a linear probabilistic model proximal splitting methods.
% Model = ml_trainproximal(Trials, Targets, Lambdas, Options...)
%
% This function allows to implement linear or logistic regression using a variety of regularization
% terms and combinations thereof using proximal splitting [1].
%
% In:
%   Trials       : training data matrix, as in ml_train
%
%   Targets      : 1d target variable vector, as in ml_train
%
%   LossType     : loss function to be used,
%                  'logistic' for classification (default)
%                  'squared' for regression
%                  'hyperbolic-secant' special-purpose for super-Gaussian estimation
%
%   Regularizers : Definition of the regularization terms. Any combination of terms is permitted.
%
%
%   Options  : optional name-value parameters to control the training details:
%
%               'regweights' : Weights of the regularizers. This is a vector of (relative) regularization
%                              parameters. If [] set to 1/N for N regularizers. (default: [])
%                              Can also be a cell array of (normalized) weight vectors; in this case 
%                              it it simultaneously optimized together with the lambdas.
%
%               'solverOptions' : cell array of name-value pairs to control how the outer ADMM solver
%                                 behaves
%
%                    'abs_tol' : Absolute tolerance criterion. (default: 10e-4)
%
%                    'rel_tol' : Relative tolerance criterion. (default: 10e-3)
%
%                    'maxit' : Maximum number of iterations. (default: 1000)
%
%                    'rho' : Initial coupling parameter. For proximal algorithms this is the coupling strength
%                            between the terms between updates. Increasing this can improve the convergence
%                            speed but too strong values can prevent any convergence. (default: 1)
%
%                    'rho_update' : Update Rho. Whether to update rho dynamically according to 3.4.1 in [2].
%                                   Note, this can sometimes cause r_norm, s_norm to "blow up". (default: true)
%
%                    'rho_cutoff' : Rho update threshold. (default: 10)
%
%                    'rho_incr' : Rho update increment factor. (default: 2)
%
%                    'rho_decr' : Rho update decrement factor. (default: 2)
%
%               'lbfgsOptions' : cell array of name-value pairs to control how the inner LFBGS solver
%                                behaves
%
%                    'm' : LBFGS history length. The number of corrections to approximate the inverse
%                          hessian matrix. (default: 6)
%
%                    'epsilon' : Tolerance criterion.  A minimization terminates when ||g|| < epsilon*max(1,||x||).
%                                (default: 1e-3)
%
%                    'past' : Distance for delta-based convergence test. (default: 0)
%
%                    'delta' : Delta for convergence test. (default: 1e-5)
%
%                    'MaxIter' : Maximum number of iterations. (default: 10)
%
%                    'linesearch' : The line search algorithm. Can be any of the following:
%                                   {'more_thuente','backtracking_armijo','backtracking_wolfe','backtracking_strong_wolfe'}
%                                   (default: more_thuente)
%
%                    'max_linesearch' : Maximum number of trials for the line search. (default: 40)
%
%                    'min_step' : Minimum step of the line search. (default: 1e-20)
%
%                    'max_step': Maximum step of the line search routine. (default: 1e20)
%
%                    'ftol' : Line search tolerance F. A parameter to control the accuracy of the
%                             line search routine. (default: 1e-4)
%
%                    'wolfe' : Coefficient for the Wolfe condition. (default: 0.9)
%
%                    'gtol' : Line search tolerance G. A parameter to control the accuracy of the
%                             line search routine. (default: 0.9)
%
%                    'xtol' : Machine precision for floating-point values. (default: 1e-16)
%
%                    'DerivativeCheck' : Derivative check using finite differences. (default: 'off')
%
%                    'Display' : Options for displaying progress. (default: 'none')
%
%               'lambdaSearch' : cell array of name-value pairs governing the regularization path search
%
%                   'lambdas' : Regulariation parameters. Controls the sparsity/simplicity of the result.
%                               Typically, this is an interval to scan. (default: 2.^(3:-0.25:-5))
%
%                   'nfolds' : Cross-validation folds. The cross-validation is used to determine the best
%                              regularization parameter value (default: 5)
%
%                   'foldmargin' : Margin (in trials) between folds. This is the number of trials omitted
%                                  between training and test sets. (default: 5)
%
%                   'cvmetric' : metric to use for parameter optimization; can be any of those supported by
%                                ml_calcloss. In particular, 'auc' is a good idea if the classification
%                                task is between highly imbalanced classes. (default: '' = auto-determine)
%
%                   'return_regpath' : Return the entire regularization path. If false, only the best model will
%                                      be returned. (default: true)
%
%
%               'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std')
%
%               'data_weights': dataset weights; optional vector of weights for each task in the
%                               training data (one element per task in a multi-task learning setting) (default: [])
%
%               'includebias': whether to include a bias param (default: true)
%
%               'verbosity': verbosity level, 0-3 (0=no output)
%
% Out:
%   Models   : a predictive model
%
% Examples:
%
% Notes:
%   When linear operators and shapes are given as string expressions the variables a to h can be used as short-hands
%   to refer to the number of array elements along respective dimension.
%
% See also:
%   ml_predictproximal
%
% References:
%  [1] Patrick L. Combettes & Jean-Christophe Pesquet, "Proximal Splitting Methods in Signal Processing",
%      in Fixed-Point Algorithms for Inverse Problems in Science and Engineering, Springer Optimization and Its Applications
%      pp. 185-212, 2011
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-02-04

% definition of regularization terms
regularizer_params = @(name) arg_subswitch({lower(name),name},{'none'},{ ...
    'none', {}, ...
    'l1', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted the norm (e.g., to impose certain types of priors).'), ...
    }, ...
    'l2', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'nonorthogonal_transform','NonorthogonalTransform'},false,[],'Linear operator is non-orthogonal. In this case an iterative method will be used that is faster and numerically more robust than letting ADMM do it.'), ...
    arg({'y','TargetValues'},[],[],'Recenter the norm around target values. This allows for regression problems as side assumptions.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    }, ...
    'l1/l2', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'g_d','GroupIndices'},uint32([]),[],'Feature group indices. This is a vector of indices that form the support of all groups; can also be a matrix. If empty, this defaults to columnwise sparsity.'), ...
    arg({'g_t','GroupSizes'},uint32([]),[],'Feature group sizes. This is a vector of successive range lengths on the indices.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    arg({'weights1','GroupWeights'},[],[],'Weights on the groups. Allows for a reweighted the norm (e.g., to impose certain types of priors).') ...
    }, ...
    'l1/linf', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'g_d','GroupIndices'},uint32([]),[],'Feature group indices. This is a vector of indices that form the support of all groups; can also be a matrix. If empty, this defaults to columnwise sparsity.'), ...
    arg({'g_t','GroupSizes'},uint32([]),[],'Feature group sizes. This is a vector of successive range lengths on the indices.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    arg({'weights1','GroupWeights'},[],[],'Weights on the groups. Allows for a reweighted the norm (e.g., to impose certain types of priors).') ...
    }, ...
    'tv2d', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.','shape','row'), ...
    arg({'useGPU','UseGPU'},false,[],'Use GPU acceleration. This is experimental and requires that UnLocBox is started with GPU support enabled.'), ...
    }, ...
    'tv3d', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.','shape','row'), ...
    arg({'useGPU','UseGPU'},false,[],'Use GPU acceleration. This is experimental and requires that UnLocBox is started with GPU support enabled.'), ...
    }, ...
    'trace', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle. When defining the linear operator as an anonymous function, the variables a to h can be used to refer to the sizes of the first 8 dimensions of x.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.','shape','row'), ...
    }}, 'Regularization term. Defines a term in the optimization problem; multiple types are supported and can be mixed freely.');

expose_handles(@solve_regularization_path,varargin{:});

arg_define([0 2],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'loss','LossType'}, 'logistic', {'logistic','squared'}, 'Loss function to be used. The logistic loss is suited for classification problems, whereas the squared loss is suited for regression problems.'), ...
    arg_sub({'regularizers','Regularizers'},{},{ ...
        regularizer_params('Term1'), ...
        regularizer_params('Term2'), ...
        regularizer_params('Term3'), ...
        regularizer_params('Term4'), ...
        regularizer_params('Term5'), ...
        regularizer_params('Term6'), ...
        regularizer_params('Term7')}, 'Definition of the regularization terms. Any combination of terms is permitted.'), ...
    arg({'regweights','TermWeights'},{[]},[],'Weights of the regularizers. This is a cell array of vectors of (relative) regularization parameters. Default is 1/N for N regularization terms. The cell array lists all possible assignments to search over.','type','expression','shape','row'), ...
    arg_sub({'solverOptions','SolverOptions'},{},{ ...
        arg({'maxit','MaxIterations'},2000,uint32([1 100 5000 10000]),'Maximum number of iterations.'), ...
        arg({'rel_tol','RelativeTolerance'},1e-3,[0 1],'Relative tolerance criterion. If the relative difference between two successive iterates is lower than this value the algorithm terminates.'),...
        arg({'abs_tol','AbsoluteTolerance'},0.000001,[0 Inf],'Absolute tolerance criterion. If the objective function value falls below this the algorithm terminates.'), ...
        arg({'rho','CouplingParameter'},4,[0 1 30 Inf],'Initial coupling parameter. For proximal algorithms this is the coupling strength between the terms between updates. Increasing this can improve the convergence speed but too strong values can prevent any convergence.'), ...
        arg({'rho_update','RhoUpdate'},true,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]. Note, this can sometimes cause r_norm, s_norm to "blow up".'), ...
        arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[0 2 20 Inf],'Rho update threshold.','guru',true), ...
        arg({'rho_incr','RhoUpdateIncr'},2.0,[1 1.5 3 Inf],'Rho update increment factor.','guru',true), ...
        arg({'rho_decr','RhoUpdateDecr'},2.0,[1 1.5 3 Inf],'Rho update decrement factor.','guru',true), ...
        arg({'warmstart','Warmstart'},true,[],'Warm-start through regularization path. Enabling this is more efficient but convergence issues can be harder to trace down.') ...
    }, 'Controls the behavior of the ADMM optimization algorithm.'), ...
    arg_sub({'lbfgsOptions','LBFGSOptions'},{},{ ...
        arg({'MaxIter','MaxIterations'},10,uint32([1 1 20 1000]),'Maximum number of iterations.'), ...
        arg({'m','HessianHistory'},6,uint32([1 4 10 20]),'LBFGS history length. The number of corrections to approximate the inverse hessian matrix.','guru',true), ...
        arg({'epsilon','Epsilon'},1e-3,[],'Tolerance criterion.  A minimization terminates when ||g|| < epsilon*max(1,||x||).'), ...
        arg({'past','DeltaDistance'},0,[],'Distance for delta-based convergence test.','guru',true), ...
        arg({'delta','Delta'},1e-5,[],'Delta for convergence test.','guru',true), ...
        arg({'linesearch','LineSearchAlgorithm'},'more_thuente',{'more_thuente','backtracking_armijo','backtracking_wolfe','backtracking_strong_wolfe'},'The line search algorithm.','guru',true), ...
        arg({'max_linesearch','MaxLineSearch'},40,uint32([0 10 100 1000]),'Maximum number of trials for the line search.','guru',true), ...
        arg({'min_step','MinStepsize'},1e-20,[],' Minimum step of the line search.','guru',true), ...
        arg({'max_step','MaxStepsize'},1e20,[],' Maximum step of the line search routine.','guru',true), ...
        arg({'ftol','FTolerance'},1e-4,[],'Line search tolerance F. A parameter to control the accuracy of the line search routine.','guru',true), ...
        arg({'wolfe','WolfeCoefficient'}, 0.9,[],'Coefficient for the Wolfe condition.','guru',true), ...
        arg({'gtol','GTolerance'},0.9,[],'Line search tolerance G. A parameter to control the accuracy of the line search routine.','guru',true), ...
        arg({'xtol','XTolerance'},1e-16,[],'Machine precision for floating-point values.','guru',true), ...
        arg({'useGPU','UseGPU'},true,[],'Run on the GPU if possible.'), ...
        arg('DerivativeCheck','off',{'off','on'},' Derivative check using finite differences.','guru',true), ...
        arg({'Display','Verbosity'},'none',{'none','on'},'Options for displaying progress.') ...
    },'Options of the inner LBFGS solver. This is for the logistic objective function.'), ...
    arg_sub({'lambdaSearch','LambdaSearch'},{},{ ...
        arg({'lambdas','Lambdas'}, 2.^(3:-0.66:-8), [0 2^-8 2^15 Inf], 'Regulariation parameters. Controls the sparsity/simplicity of the result. Typically, this is an interval to scan, such as 2.^(10:-1:-15).'), ...
        arg({'nfolds','NumFolds'},5,[],'Cross-validation folds. The cross-validation is used to determine the best regularization parameter. If in 0..1, k calulated as fraction of #trials, if in -1..0, taken as the p in p-holdout, if given as [low, high], taken as a fractional interval for interval holdout.','shape','row'),...
        arg({'force_cv','ForceCV'},false,[],'Force cross-validation. This will perform a nested cross-validation even if there is only one lambda/regweight (i.e., a search wouldn''t be strictly necessary). Can be useful to simultaneously run multiple within-task/subject cross-validations given fixed reg parameters.','shape','row'),...
        arg({'foldmargin','FoldMargin'},0,uint32([0 0 10 1000]),'Margin between folds. This is the number of trials omitted between training and test set.'), ...
        arg({'cvmetric','ParameterMetric'},'',{'','kld','nll','mcr','mae','mse','max','rms','bias','medse','auc','cond_entropy','cross_entropy','f_measure'},'Metric for Parameter Optimization. By default auto-determined; can be any of the ml_calcloss-supported metrics. In particular, auc is a good idea if the classification task is between highly imbalanced classes.') ...
        arg({'return_regpath','ReturnRegpath'}, true, [], 'Return the entire regularization path. This is for the best relative weighting of terms. If false, only the best model will be returned.'), ...
        arg({'return_reggrid','ReturnReggrid'}, false, [], 'Return the entire regularization grid. This also returns regularization paths for all other relative weightings. Warning: this can require a lot of memory (depending on model size).'), ...
        arg({'history_traces','HistoryTraces'}, false, [], 'Return history traces. If true, optimization history traces will be returned. Warning: this will require a very large amount of memory (depending on model size).'), ...
    }, 'Controls the search for the optimal regularization parameter.'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg_nogui({'shape','Shape'}, [], [], 'Reshaping for features. Allows to reshape (perhaps vectorized) features into a particular representation.','shape','row'), ...
    arg({'data_weights','DataWeights'}, [], [], 'Dataset weights. Optional vector of weights for each task in the training data (one element per task in a multi-task learning setting).'), ...
    arg({'verbosity','Verbosity'},1,uint32([1 5]),'Diagnostic output level. Zero is off, 1 only shows cross-validation diagnostics, 2 shows solver diagnostics, 3 shows iteration diagnostics.'), ...
    arg({'continuous_targets','ContinuousTargets','Regression'}, false, [], 'Whether to use continuous targets. This allows to implement some kind of damped regression approach.'),...
    arg({'votingScheme','VotingScheme'},'1v1',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'), ...
    arg({'parallel_scope','ParallelScope'},[],[],'Optional parallel scope. If this is a cell array of name-value pairs, cluster resources will be acquired with these options for the duration of bci_train (and released thereafter) Options as in env_acquire_cluster.','type','expression'), ...
    arg({'engine_cv','ParallelEngine','engine'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine to use. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'includebias','IncludeBias','bias'},true,[],'Include bias param. Also learns an unregularized bias param (strongly recommended for typical classification problems).'));

if ~iscell(targets)
    trials = {trials};
    targets = {targets}; 
end

for t=1:length(trials)
    trials{t} = real(trials{t}); end

% find all target classes (if classification)
nTasks = length(targets);
classes = unique(vertcat(targets{:}));
if length(classes) > 2 && strcmp(loss,'logistic') && ~continuous_targets
    % in the multi-class case we use the voter for now (TODO: use softmax loss instead)
    model = ml_trainvote(trials, targets, votingScheme, @ml_trainproximal, @ml_predictproximal, varargin{:});
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
        
    if isscalar(lambdaSearch.nfolds)
        % if nfolds is in [0..1], we take it as a function of #trials
    if lambdaSearch.nfolds < 1 && lambdaSearch.nfolds > 0
        lambdaSearch.nfolds = round(lambdaSearch.nfolds*mean(cellfun('length',targets))); end
        % (if instead it's in [-1..0], we take it as the p in p-holdout)
        % (else if an integer, we take it as the number of folds)
    nFolds = ceil(abs(lambdaSearch.nfolds));
    elseif isequal(size(lambdaSearch.nfolds), [1 2])
        % if nfolds is of the form [low high], we treat these numbers as an interval specification,
        % where low and high are taken as fractions of the number of trials
        nFolds = 1;
    else
        error('NumFolds format is unsupported.');
    end
    
    % sanitize some more inputs
    solverOptions.verbose = max(0,verbosity-1);    
    if isnumeric(regweights)
        regweights = {regweights}; end
    nRegweights = length(regweights);
    
    % lambdas need to be sorted in descending order for the warm-starting to work
    nLambdas = length(lambdaSearch.lambdas);
    lambdaSearch.lambdas = sort(lambdaSearch.lambdas,'ascend');
    if strcmp(lambdaSearch.cvmetric,'mcr')
        lambdaSearch.cvmetric = ''; end
    
    % determine featureshape and vectorize data if necessary 
    [featureshape,trials,vectorize_trials] = utl_determine_featureshape(trials,shape);
    weightshape = [featureshape nTasks];
    
    % optionally scale the data
    sc_info = hlp_findscaling(vertcat(trials{:}),scaling);
    trials = cellfun(@(t)hlp_applyscaling(t,hlp_findscaling(t,scaling)),trials,'UniformOutput',false);
    
    % optionally remap target labels to -1,+1
    if strcmp(loss,'logistic') && length(classes) == 2 && ~continuous_targets
        for t=1:nTasks
            targets{t}(targets{t}==classes(1)) = -1;
            targets{t}(targets{t}==classes(2)) = +1;
        end
    end
    
    % ensure that data_weights exists and is scaled properly (we normalize data_weights to sum to
    % nTasks, since the regularizers will usually also be scaled by nTasks)
    if isempty(data_weights)
        data_weights = ones(1,nTasks); end
    data_weights = data_weights/sum(data_weights)*nTasks;
        
    % learn a sequence of models across the given lambda's, on all the data (i.e. the regularization path)
    if verbosity
        disp('Running optimization...'); end
    
    
    % run a cross-validation to score the lambdas and regweights

    % loss_means{regweight,task}(fold,lambda) is the average loss for a given task, regularization weight setting, cross-validation fold, and lambda setting
    loss_means = cell(nRegweights,nTasks);    
    % predictions{regweight,task}(trial,lambda) is the classifier prediction for a given task, regweight setting, trial, and lambda choice
    predictions = repmat(cellfun(@(t)zeros(length(t),nLambdas),targets(:)','UniformOutput',false),nRegweights,1);
    % foldid{task}(trial) is the fold in which a given trial is in the test set, for a given task
    if isequal(size(lambdaSearch.nfolds), [1 2])
        % interval form
        p = lambdaSearch.nfolds;
        foldids = cellfun(@(t)(0:length(t)-1)/length(t) > p(1) & (0:length(t)-1)/length(t) < p(2),targets,'UniformOutput',false);
    elseif lambdaSearch.nfolds < 0 && lambdaSearch.nfolds > -1
        p = abs(lambdaSearch.nfolds);
        % negative fractional value encodes p-holdout (positive fractional value is already defined
        % as the a fraction of the number of trials)
        foldids = cellfun(@(t)(0:length(t)-1)/length(t)>(1-p),targets,'UniformOutput',false);
    else
        foldids = cellfun(@(t)1+floor((0:length(t)-1)/length(t)*nFolds),targets,'UniformOutput',false);
    end
    
    if (nLambdas*nRegweights) > 1 || lambdaSearch.force_cv
        % for each fold...
        model_seq = cell(nFolds,nRegweights); % model_seq(fold,regweight}{lambda}{task} is the model for a given fold, regweight and lambda setting, and task
        history_seq = cell(nFolds,nRegweights); % history_seq(fold,regweight}{lambda}( is a struct of optimization histories for a given fold, regweight and lambda setting, for all concurrent tasks
        jobs = {}; % compute jobs
        for f = 1:nFolds
            % determine training and test set masks
            % TODO: calc all this per fold and don't recalc below
            testmask{f} = cellfun(@(foldid)foldid==f,foldids,'UniformOutput',false); % testmask{fold}{task}(trial) a bitmask of test-set trials for a given task
            trainmask{f} = cellfun(@(x)~x,testmask{f},'UniformOutput',false);           % trainmask{fold}{task}(trial) is a bitmask of train-set trials
            % cut train/test margins into trainmask
            for t=1:nTasks
                testpos = find(testmask{f}{t});
                for j=1:lambdaSearch.foldmargin
                    trainmask{f}{t}(max(1,testpos-j)) = false;
                    trainmask{f}{t}(min(length(testmask{f}{t}),testpos+j)) = false;
                end
            end

            % set up design matrices
            [A{f},y{f},B{f},z{f}] = deal(cell(1,nTasks));
            for t=1:nTasks
                % training data
                A{f}{t} = trials{t}(trainmask{f}{t},:);
                y{f}{t} = targets{t}(trainmask{f}{t});
                % test data
                B{f}{t} = [trials{t}(testmask{f}{t},:) ones(nnz(testmask{f}{t}),double(includebias))];
                z{f}{t} = targets{t}(testmask{f}{t});
            end

            % for each relative regularization term weighting...
            for w = 1:nRegweights
                jobs{end+1} = {@hlp_diskcache,'predictivemodels',@solve_regularization_path,A{f},y{f},lambdaSearch.lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,regweights{w},weightshape,data_weights}; end            
        end
        
        % run the jobs
        results = par_schedule(jobs, 'engine',engine_cv, 'scope',parallel_scope);

        % evaluate results
        predictions = repmat(cellfun(@(t)zeros(length(t),nLambdas),targets(:)','UniformOutput',false),nRegweights,1);
        ji = 1;
        for f = 1:nFolds
            for w = 1:nRegweights
                [model_seq{f,w},history_seq{f,w}] = deal(results{ji}.regpath,results{ji}.hist); ji = ji+1;
                % for each task...
                for t = 1:nTasks
                    % calc test-set predictions for each model
                    for m=nLambdas:-1:1
                        predictions{w,t}(testmask{f}{t},m) = (B{f}{t}*model_seq{f,w}{m}{t}(:))'; end
                    if strcmp(loss,'logistic')
                        predictions{w,t}(testmask{f}{t},:) = 2*(1 ./ (1 + exp(-predictions{w,t}(testmask{f}{t},:))))-1; end

                    % evaluate test-set losses
                    if isempty(lambdaSearch.cvmetric)
                        if strcmp(loss,'logistic')
                            loss_means{w,t}(f,:) = mean(~bsxfun(@eq,z{f}{t},sign(predictions{w,t}(testmask{f}{t},:))));
                        else
                            loss_means{w,t}(f,:) = mean((bsxfun(@minus,z{f}{t},predictions{w,t}(testmask{f}{t},:))).^2);
                        end
                    else
                        for m=1:nLambdas
                            loss_means{w,t}(f,m) = ml_calcloss(lambdaSearch.cvmetric,z{f}{t},predictions{w,t}(testmask{f}{t},m)); end
                    end
                end
            end
        end
    else
        % we skip the nested cross-validation if there is only one lambda and one regweight
        disp('Skipping nested cross-validation (only 1 lambda/regweight)...'); 
        model_seq = cell(nFolds,nRegweights); 
        history_seq = cell(nFolds,nRegweights);
        loss_means = repmat({zeros(nFolds,nLambdas)},[nRegweights,nTasks]);
        lambdaSearch.return_regpath = false;
    end
    
    % pick best lambda and regweights across tasks (averaging over folds)
    losses = zeros(nRegweights,nLambdas,nTasks);
    for t=1:nTasks
        for w=1:nRegweights
            losses(w,:,t) = mean(loss_means{w,t}, 1); end        
        joint_losses = mean(losses,3); % nRegweights x nLambdas
        % find per-task best lambda/regweights
        [best_regweight_indices,best_lambda_indices] = find(losses(:,:,t) == min(vec(losses(:,:,t))));
        [dummy,idx] = max(best_lambda_indices); %#ok<ASGLU>
        best_regweights{t} = regweights{best_regweight_indices(idx)}(:)';
        best_lambda{t} = lambdaSearch.lambdas(best_lambda_indices(idx));
    end
    % find all optima in the loss surface and then pick the minimum at highest lambda (if multiple)
    [best_regweight_indices,best_lambda_indices] = find(joint_losses == min(joint_losses(:)));
    [dummy,idx] = max(best_lambda_indices); %#ok<ASGLU>
    joint_best_regweights = regweights{best_regweight_indices(idx)}(:)';
    joint_best_lambda = lambdaSearch.lambdas(best_lambda_indices(idx));

    % pick the model at the minimum...
    if lambdaSearch.return_regpath
        % run the whole regularization path for the jointly best regweight combination
        res = hlp_diskcache('predictivemodels',@solve_regularization_path,trials,targets,lambdaSearch.lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,joint_best_regweights,weightshape,data_weights);
        [regpath,history] = deal(res.regpath, res.hist);
        model.regularization_path = regpath;                                  % the model for a given {lambda}{task} at best regweights, for whole data
        model.regularization_loss = permute(losses(best_regweight_indices(idx),:,:),[2,3,1]);   % the associated loss estimatses for a given (lambda,task)
        model.ws = regpath{find(lambdaSearch.lambdas == joint_best_lambda,1)};% the best model for a given {task}
    else
        res = hlp_diskcache('predictivemodels',@solve_regularization_path,trials,targets,joint_best_lambda,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,joint_best_regweights,weightshape,data_weights);
        [tmp,history] = deal(res.regpath, res.hist);
        model.ws = tmp{1};  % optimal model for each {task}
    end
    
    model.w = model.ws;
    if length(model.w) == 1
        model.w = model.w{1}; end                   % optimal model for each {task}, or the model if only one task given (without the enclosing cell array)
    
    if lambdaSearch.return_reggrid
        model.regularization_grid = model_seq; end  % sequence of models for each {fold,regweight}{lambda}{task}
    if lambdaSearch.history_traces
        model.fold_history = history_seq;           % structure of regpath history for each {fold,regweight}{lambda} -- HUGE!
        model.regularization_history = history;     % structure of regpath history at best lambda/regweights
    end
    model.loss_means = loss_means;                  % overall loss for each {regweight,task}(fold,lambda)
    model.task_losses = losses;                     % average loss for each (reweight,task,lambda)
    model.joint_losses = joint_losses;              % average loss for each (regweight,lambda)
    model.best_regweights = best_regweights;        % best regweights for each {task}
    model.best_lambda = best_lambda;                % best lambda for each {task}
    model.joint_best_regweights = joint_best_regweights; % best regweights from all tasks
    model.joint_best_lambda = joint_best_lambda;         % best lambda from all tasks
    model.classes = classes;                        % set of class labels in training data
    model.continuous_targets = continuous_targets;  
    model.includebias = includebias;                % whether a bias is included in the model
    model.vectorize_trials = vectorize_trials;      % whether trials need to be vectorized first
    model.featureshape = featureshape;              % shape vector for features
    model.sc_info = sc_info;                        % overall scaling info
    model.loss = loss;                              % loss function name
end



% learn the regularization path
function res = solve_regularization_path(A,y,lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizersArg,regweights,weightshape,data_weights)
% solve_regularization_path_version<1.0.3>
if ~includebias
    error('This implementation currently requires that a bias is included.'); end
nTasks = length(A);

% m trials, n features, per task
m = cellfun('size',A,1);
n = cellfun('size',A,2);
if length(unique(n)) > 1
    error('Each task must have the same number of features.'); end

% w is the concatenation of model weights for all tasks, followed by the unregularized biases for each task
w = zeros(sum(n) + nTasks,1);

% set up the design matrix A & label vector y
A = cellfun(@(A)double(A),A,'UniformOutput',false);
y = cellfun(@(y)double(y(:)),y,'UniformOutput',false);

% set up the data-dependent loss function to use
switch loss
    case 'logistic'
        C = cellfun(@(A,y)[bsxfun(@times,-y,A) -y],A,y,'UniformOutput',false);
        if lbfgsOptions.useGPU
            try
                C = cellfun(@gpuArray,C,'UniformOutput',false);
            catch e
                disp_once(['Could not enable GPU support: ' e.message]);
            end
        end
        Cp = cellfun(@transpose,C,'UniformOutput',false);
        lossfunc.prox = @(x,gamma,x0) prox_logistic_multitask(C,Cp,x,gamma,x0,m,n,hlp_struct2varargin(lbfgsOptions),data_weights);
        lossfunc.eval = @(x,lambda) lambda*obj_logistic_multitask(C,x,m,n,data_weights);
    case 'squared'
        % append a bias to the design matrix
        if length(A)>1
            error('Squared loss with for multi-task case not yet fully implemented.'); end
        Ao = cellfun(@(A)[A ones(size(A,1),1)],A,'UniformOutput',false);
        mm = cellfun('size',Ao,1);
        nn = cellfun('size',Ao,2);
        % choose the right prox operator
        if solverOptions.rho_update
            lossfunc.prox = @(x,gamma,x0) prox_squared_iter(Ao{1},y{1},1/gamma,x,zeros(size(x)),mm,nn,x0);
        else            
            Atb = cellfun(@(Ao,y)Ao'*y,Ao,y,'UniformOutput',false);
            [L,U] = deal(cell(1,nTasks));
            for t=1:nTasks
                [L{t},U{t}] = factor(Ao{t},solverOptions.rho/data_weights(t)); end
            lossfunc.prox = @(x,gamma,x0) prox_squared_factored_multitask(Ao,Atb,L,U,solverOptions.rho,x,zeros(size(x)),mm,nn,data_weights);
        end
        lossfunc.eval = @(x,lambda) lambda*obj_squared_multitask(Ao,y,x,mm,nn,data_weights);
    case 'hyperbolic-secant'
        % lossfunc.prox = @(x,gamma,x0) prox_hs(y,x,gamma,x0,hlp_struct2varargin(lbfgsOptions));
        % lossfunc.eval = @(x,lambda) lambda*obj_hs(x,b);
        error('Hyperbolic-secant loss is not yet implemented.');
    otherwise
        error('Unsupported loss function.');
end
lossfunc.y0 = [];
lossfunc = @(lambda)setfield(setfield(lossfunc,'prox',@(x,gamma,x0)lossfunc.prox(x,gamma*lambda,x0)),'eval',@(x)lossfunc.eval(x,lambda)); %#ok<SFLD>


% ensure that regularizers is a cell array of structs
regularizers = {};
if isstruct(regularizersArg)
    for k=1:length(fieldnames(regularizersArg))
        if isfield(regularizersArg,['term' num2str(k)])
            regularizers{end+1} = regularizersArg.(['term' num2str(k)]); end %#ok<AGROW>
    end
else
    regularizers = regularizersArg;
end


% set up the regularization functions one by one
regfuncs = {};
for t = 1:length(regularizers)
    param = regularizers{t};
    if ~strcmp(param.arg_selection,'none')
        regfunc = struct();
        if isfield(param,'weights')
            param.weights2 = param.weights; end
        param.verbose = max(0,solverOptions.verbose-2);
        
        % rename & evaluate the linear operator expressions
        if ischar(param.A)
            try
                [a,b,c,d,e,f,g,h] = size(reshape(w(1:sum(n)),weightshape)); %#ok<ASGLU>
                param.A = eval(param.A);
            catch e
                env_handleerror(e);
                disp(['This param does not evaluate correctly: '  param.A]);
            end
        end
        
        % if the linear operator happens to accept weights in the shape of the original features
        % (and the numels are matching) then we reshape the weights to that shape before applying
        % the linear operator
        try
            rA = param.A;
            rA(reshape(w(1:sum(n)),weightshape));
            param.A = @(x)rA(reshape(x(1:sum(n)),weightshape));
        catch
            try
                param.A(w(1:sum(n)));
                param.A = @(x)rA(x(1:sum(n)));
            catch e
                % sanity check: if this happens either your linear operator is incorrect or
                % you need to specify NumberOfElements for this term
                error(['The linear operator ' char(param.A) ' is not applicable to the weights w. Check for syntax errors and sizes.']);
            end
        end
        shape_A_out = size(param.A(w));
        
        % if shape for the term is unspecified we assume that it is the output shape of the A
        % operator
        if isfield(param,'shape') && isempty(param.shape)
            param.shape = shape_A_out; end
        if isfield(param,'shape') && ischar(param.shape)
            [a,b,c,d,e,f,g,h] = size(reshape(w(1:sum(n)),weightshape)); %#ok<ASGLU>
            param.shape = eval(param.shape);
        end
        
        % set the A matrix for future reference
        if isfield(param,'A')
            rA = param.A;
        else
            rA = @(x)x(1:sum(n));
        end
        
        % move the A parameter into regfunc.L (handled by ADMM)
        if isfield(param,'A') && ~(strcmp(param.arg_selection,'l2') && param.nonorthogonal_transform)
            regfunc.L = param.A;
            % remove fields from param
            param = rmfield(param,'A');
            if isfield(param,'At')
                param = rmfield(param,'At'); end
        else
            regfunc.L = @(x)x(1:sum(n));
        end
        
        % now turn .L into a matrix (since we actually need it in matrix form)
        % the calculation is cached since it's quite slow for large parameter spaces
        regfunc.L = operator_to_matrix(regfunc.L,numel(w));
        
        regfunc.param = param;
        
        vec = @(x)x(:);
        if isfield(param,'weights')
            if isempty(param.weights)
                param.weights = ones(prod(shape_A_out),1); end
            shaped_weights = reshape(param.weights,shape_A_out);
        else
            shaped_weights = ones(shape_A_out);
        end
        switch param.arg_selection
            case 'l1'
                if (isempty(param.weights) || all(param.weights(:)==1)) && ~isfield(param,'A')
                    regfunc.prox = @(x,gamma,x0) prox_l1_simple(x,gamma);
                    regfunc.eval = @(x,lambda) lambda*sum(abs(x));
                else
                    regfunc.prox = @(x,gamma,x0) prox_l1(x,gamma,param);
                    regfunc.eval = @(x,lambda) lambda*norm(vec(shaped_weights.*rA(x)),1);
                end
            case 'l2'
                if isempty(param.y)
                    param.y = zeros(prod(shape_A_out),1); end
                if isfield(param,'A')
                    A = operator_to_matrix(param.A,sum(n));
                    regfunc.prox = @(x,gamma,x0) prox_squared_iter(A,param.y,1/gamma,x,zeros(size(x)),sum(n),x0);
                    regfunc.eval = @(x,lambda) lambda*obj_squared(A,param.y,x(1:sum(n)));
                else                
                    regfunc.prox = @(x,gamma,x0) prox_l2_simple(x,gamma);
                    regfunc.eval = @(x,lambda) lambda*norm(shaped_weights(:).*(vec(rA(x)) - param.y(:)),2).^2;
                end
            case 'l1/l2'
                if isfield(param,'A')
                    error('The linear operator for the group sparsity prox operator is not implemented. You can however apply it by using the SDMM algorithm.'); end
                if isempty(param.g_d) && isempty(param.g_t) && (isempty(param.weights2)||all(param.weights2(:)==1)) && (isempty(param.weights1)||all(param.weights1(:)==1))
                    regfunc.prox = @(x,gamma,x0) prox_l12_simple(x,gamma,shape_A_out);
                    regfunc.eval = @(x,lambda)lambda*norm_l12_simple(rA(x),shape_A_out);
                else                    
                    if isempty(param.g_d) && isempty(param.g_t)
                        param.g_d = (1:prod(shape_A_out));
                        param.g_t = shape_A_out(1)*ones(1,prod(shape_A_out(2:end)));
                    end
                    if isempty(param.weights1)
                        param.weights1 = ones(numel(param.g_t),1); end
                    if isempty(param.weights2)
                        param.weights2 = ones(prod(shape_A_out),1); end
                    regfunc.prox = @(x,gamma,x0) prox_l12(x,gamma,param);
                    regfunc.eval = @(x,lambda) lambda*norm_l12(rA(x),param.g_d,param.g_t,param.weights2,param.weights1);
                end
            case 'l1/linf'
                if isfield(param,'A')
                    error('The linear operator for the group sparsity prox operator is not implemented. You can however apply it by using the SDMM algorithm.'); end
                if isempty(param.g_d) && isempty(param.g_t) && (isempty(param.weights2)||all(param.weights2(:)==1)) && (isempty(param.weights1)||all(param.weights1(:)==1))
                    regfunc.prox = @(x,gamma,x0) prox_l1inf_simple(x,gamma,shape_A_out);
                    regfunc.eval = @(x,lambda)lambda*norm_l1inf_simple(rA(x),shape_A_out);
                else                    
                    if isempty(param.g_d) && isempty(param.g_t)
                        param.g_d = (1:prod(shape_A_out));
                        param.g_t = shape_A_out(1)*ones(1,prod(shape_A_out(2:end)));
                    end
                    if isempty(param.weights1)
                        param.weights1 = ones(numel(param.g_t),1); end
                    if isempty(param.weights2)
                        param.weights2 = ones(prod(shape_A_out),1); end
                    regfunc.prox = @(x,gamma,x0) prox_l1inf(x,gamma,param);
                    regfunc.eval = @(x,lambda) lambda*norm_l1inf(rA(x),param.g_d,param.g_t,param.weights2,param.weights1);
                end
            case 'tv2d'
                if isfield(param,'A')
                    error('The linear operator for the total-variation prox operator is not implemented. You can however apply it by using the SDMM algorithm.'); end
                regfunc.prox = @(x,gamma,x0) prox_tv(x,gamma,param);
                regfunc.eval = @(x,lambda) lambda*tv_norm(rA(x),param.shape);
            case 'tv3d'
                if isfield(param,'A')
                    error('The linear operator for the total-variation prox operator is not implemented. You can however apply it by using the SDMM algorithm.'); end
                regfunc.prox = @(x,gamma,x0) prox_tv3d(x,gamma,param);
                regfunc.eval = @(x,lambda) lambda*tv_norm3d(rA(x),param.shape);
            case 'trace'
                regfunc.prox = @(x,gamma,x0) prox_nuclear_simple(x,gamma,shape_A_out);
                regfunc.eval = @(x,lambda) lambda*norm_nuclear_simple(rA(x),shape_A_out);
            otherwise
                error('Unrecognized regularization type requested.');
        end
        regfunc.y0 = [];
        regfuncs{end+1} = @(lambda) setfield(setfield(regfunc,'prox',@(x,gamma,x0)regfunc.prox(x,gamma*lambda,x0)),'eval',@(x)regfunc.eval(x,lambda)); %#ok<AGROW,SFLD>
    end
end

if iscell(regweights) && numel(regweights) == 1
    regweights = regweights{1}; end
if isempty(regweights)
    regweights = ones(1,length(regfuncs)); end
regweights = regweights/sum(regweights); 

% learn the regularization path
if verbosity
    disp('solving regularization path...'); end
y0 = {};
nLambdas = length(lambdas);
regpath = cell(nLambdas,1);
hist = cell(1,nLambdas);
for k =1:nLambdas
    
    % set up parameters
    termweights = [1,lambdas(k)*regweights];
    lossfunc = lossfunc(1);
    lossfunc.L = [];
    if ~isempty(y0)
        lossfunc.y0 = y0{1}; end
    lossfunc.param = struct();
    for r = 1:length(regfuncs)
        tmpregfuncs(r) = regfuncs{r}(termweights(1+r)); %#ok<AGROW>
        if ~isempty(y0)
            tmpregfuncs(r).y0 = y0{1+r}; end %#ok<AGROW>
    end
    
    % we stash the termweights in the solverOptions because the hlp_diskcache below will by default
    % not parse the tmpregfuncs anonymous function deep enough to discover the termweights, and thus 
    % cause cache collisions (this can be resolved by setting the serialize_anonymous_fully option
    % to true, but since that is a global option it could cause unexpected behavior in the rest of
    % BCILAB)
    solverOptions.termweights = termweights;
    
    % solve
    t0 = tic;
    if verbosity
        fprintf('  scanning lambda = %f (%i/%i)...',lambdas(k),k,nLambdas); end
    if solverOptions.warmstart
        [w,y0,rho,hist{k}] = hlp_diskcache('intermediate',@consensus_admm,w,[lossfunc tmpregfuncs],solverOptions); %#ok<ASGLU>
    else
        [w,y0dummy,rho,hist{k}] = hlp_diskcache('intermediate',@consensus_admm,zeros(size(w)),[lossfunc tmpregfuncs],solverOptions); %#ok<ASGLU>
    end
    if verbosity
        fprintf(' %i iters; t = %.1fs\n',length(hist{k}.objval),toc(t0)); end
    
    % assemble output weights for each task
    if includebias && nTasks > 1
        offsets = cumsum([1 n(1:end-1)]);
        for t=1:nTasks
            regpath{k}{t} = w([offsets(t)+(0:n(t)-1) end-length(m)+t]); end
    else
        regpath{k}{1} = w;
    end
end

[res.regpath,res.hist] = deal(regpath,hist);


% --- multi-task logistic loss code ---

function [val,grad] = obj_proxlogistic_multitask(x,C,Cp,z,gamma,m,n,data_weights)
% objective function for the multi-task logistic loss proximity operator (effectively l2-regularized logreg)
offsets = cumsum([1 n(1:end-1)]);
for t=length(m):-1:1
    % move bias from end to inline
    xt = x([offsets(t)+(0:n(t)-1) end-length(m)+t]);
    zt = z([offsets(t)+(0:n(t)-1) end-length(m)+t]);
    ecx = exp(C{t}*xt);
    scaling = (gamma*data_weights(t)/m(t));
    val{t} = (1/2)*sum((xt-zt).^2) + scaling*gather(sum(log1p(ecx)));
    if ~isfinite(val{t})
        ecx(~isfinite(ecx(:))) = 2.^50;
        val{t} = (1/2)*sum((xt-zt).^2) + scaling*gather(sum(log1p(ecx)));
    end
    grad{t} = (xt - zt) + scaling*gather(Cp{t}*(ecx./(1+ecx)));    
end
val = sum([val{:}]);
grad = vertcat(grad{:});
% move biases back to end
grad = [grad;grad(cumsum(n+1))]; grad(cumsum(n+1)) = [];

function x = prox_logistic_multitask(C,Cp,z,gamma,x0,m,n,args,data_weights)
x = liblbfgs(@(w)obj_proxlogistic_multitask(w,C,Cp,z,gamma,m,n,data_weights),x0,args{:});

function obj = obj_logistic_multitask(C,x,m,n,data_weights)
obj = 0;
offsets = cumsum([1 n(1:end-1)]);
for t=1:length(m)
    xt = x([offsets(t)+(0:n(t)-1) end-length(m)+t]);
    obj = obj + gather(sum(log1p(exp(C{t}*xt))))*(data_weights(t)/m(t)); 
end


% --- multi-task square loss code  ---

function x = prox_squared_factored_multitask(A,Atb,L,U,rho,z,u,m,n,data_weights)
% this version can only be used if rho stays constant (TODO: confirm the use of data_weights as correct)
scaling = rho/data_weights;
offsets = cumsum([1 n(1:end-1)]);
for t=length(m):-1:1
    q = Atb{t} + scaling(t)*(z([offsets(t)+(0:n(t)-1) end-length(m)+t]) - u([offsets(t)+(0:n(t)-1) end-length(m)+t]));
    if(m(t) >= n(t))
        x{t} = U{t}\(L{t}\q);
    else
        x{t} = q/scaling(t) - (A{t}'*(U{t}\(L{t}\(A{t}*q))))/scaling(t)^2;
    end
end
x = vertcat(x{:});
x = [x;x(cumsum(n+1))]; x(cumsum(n+1)) = [];
    
function obj = obj_squared_multitask(A,b,x,m,n,data_weights)
obj = 0;
for t=length(m):-1:1
    xt = x([offsets(t)+(0:n(t)-1) end-length(m)+t]); 
    obj = obj + data_weights(t)*0.5*norm(A*xt - b,2).^2;
end


% --- logistic loss code ---

function [val,grad] = obj_proxlogistic(x,C,Cp,z,gamma,m)
% objective function for the logistic loss proximity operator (effectively l2-regularized logreg)
ecx = exp(C*x);
scaling = (gamma/m);
val = (1/2)*sum((x-z).^2) + scaling*gather(sum(log1p(ecx)));
if ~isfinite(val)
    ecx(~isfinite(ecx(:))) = 2.^50;
    val = (1/2)*sum((x-z).^2) + scaling*gather(sum(log1p(ecx)));
end
grad = (x - z) + scaling*gather(Cp*(ecx./(1+ecx)));

function x = prox_logistic(C,Cp,z,gamma,x0,m,args)
x = liblbfgs(@(w)obj_proxlogistic(w,C,Cp,z,gamma,m),x0,args{:});

function obj = obj_logistic(C,x,m)
obj = gather(sum(log1p(exp(C*x))))/m;

           
% --- square loss code  ---

function x = prox_squared_factored(A,Atb,L,U,rho,z,u,m,n)
% this version can only be used if rho stays constant
q = Atb + rho*(z - u);
if(m >= n)
    x = U\(L\q);
else
    x = q/rho - (A'*(U\(L\(A*q))))/rho^2;
end

function x = prox_squared_iter(A,b,rho,z,u,n,x0)
[x, flag, relres, iters] = lsqr([A; sqrt(rho)*speye(n)], [b; sqrt(rho)*(z-u)], [], [], [], [], x0); %#ok<NASGU,ASGLU>

function obj = obj_squared(A,b,x)
obj = 0.5*norm(A*x - b,2).^2;


% --- some useful prox operators & norms ---

function x = prox_l1_simple(z, gamma)
% for the l1 norm
x = max(0,z-gamma) - max(0,-z-gamma);

function x = prox_l12_simple(z, gamma, shape)
% for the columnwise group l1/l2 norm
z = reshape(z,shape);
x = bsxfun(@times,max(0,1-gamma./sqrt(sum(z.^2))),z);
x = x(:);

function x = prox_l2_simple(z, gamma)
% for the l2 norm
x = bsxfun(@times,max(0,1-gamma./sqrt(sum(z.^2))),z);

function x = prox_l1inf_simple(z, gamma, shape)
% for the columnwise group l1/linf norm
z = reshape(z,shape);
x = bsxfun(@times,max(0,1-gamma/max(abs(z))),z);
x = x(:);

function x = prox_nuclear_simple(z, gamma, shape)
% for the trace norm on first 2 dimensions
z = reshape(z,shape);
if ndims(z)>2
    siz = size(z);
    z = reshape(z,siz(1),siz(2),[]);
    for k=1:size(z,3)
        [U,S,V] = svd(z(:,:,k),'econ');
        S = diag(max(0,diag(S)-gamma));
        z(:,:,k) = U*S*V.';
    end
    x = reshape(z,siz);
else
    [U,S,V] = svd(z,'econ');
    S = diag(max(0,diag(S)-gamma));
    x = U*S*V.';
end
x = x(:);

function n = norm_l12_simple(z, shape)
% for the columnwise group l1/l2 norm
n = sqrt(sum(reshape(z.^2,shape)));
n = sum(n(:));

function n = norm_l1inf_simple(z, shape)
% for the columnwise group l1/linf norm
n = max(reshape(abs(z),shape));
n = sum(n(:));

function n = norm_nuclear_simple(z, shape)
% for the trace norm on first 2 dimensions
z = reshape(z,shape);
if ndims(z)>2
    siz = size(z);
    z = reshape(z,siz(1),siz(2),[]);
    n = 0;
    for k=1:size(z,3)
        n = n+sum(svd(z(:,:,k))); end
else
    n = sum(svd(z));
end


% -- hyperbolic secant distribution loss code ---

function [val,grad] = obj_proxhs(x,b,z,gamma)
% objective function for the hyperbolic-secant loss proximity operator
zz = x-b;
mz = abs(zz);
ezzmz = exp(zz-mz);
enzzmz = exp(-zz-mz);
val = (1/2)*norm(x - z).^2 + gamma*sum(mz + log(ezzmz+enzzmz)-log(2/pi));
grad = (x - z) + gamma*(ezzmz-enzzmz)./(ezzmz+enzzmz);

function x = prox_hs(b,z,gamma,x0,args)
x = liblbfgs(@(w)obj_proxhs(w,b,z,gamma),x0,args{:});

function obj = obj_hs(x,b)
x = x-b;
ax = abs(x);
obj = sum(mz + log(exp(x-ax)+exp(-x-ax))-log(2/pi));


% --- helper functions ---

function [L,U] = factor(A, rho)
% note: rho is 1/gamma
[m,n] = size(A);
if (m >= n)
    L = chol(A'*A + rho*speye(n),'lower');
else
    L = chol(speye(m) + 1/rho*(A*A'),'lower');
end
% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');


