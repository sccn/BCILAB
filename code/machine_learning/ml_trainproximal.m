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
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted the norm (e.g., to impose certain types of priors).'), ...
    }, ...
    'l2', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'nonorthogonal_transform','NonorthogonalTransform'},false,[],'Linear operator is non-orthogonal. In this case an iterative method will be used that is faster and numerically more robust than letting ADMM do it.'), ...
    arg({'y','TargetValues'},[],[],'Recenter the norm around target values. This allows for regression problems as side assumptions.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    }, ...
    'l1/l2', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'g_d','GroupIndices'},[],[],'Feature group indices. This is a vector of indices that form the support of all groups; can also be a matrix. If empty, this defaults to columnwise sparsity.'), ...
    arg({'g_t','GroupSizes'},[],[],'Feature group sizes. This is a vector of successive range lengths on the indices.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    arg({'weights1','GroupWeights'},[],[],'Weights on the groups. Allows for a reweighted the norm (e.g., to impose certain types of priors).') ...
    }, ...
    'l1/linf', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'g_d','GroupIndices'},[],[],'Feature group indices. This is a vector of indices that form the support of all groups; can also be a matrix. If empty, this defaults to columnwise sparsity.'), ...
    arg({'g_t','GroupSizes'},[],[],'Feature group sizes. This is a vector of successive range lengths on the indices.'), ...
    arg({'weights','FeatureWeights'},[],[],'Weights on the features. Allows for a reweighted norm (e.g., to impose certain types of priors).'), ...
    arg({'weights1','GroupWeights'},[],[],'Weights on the groups. Allows for a reweighted the norm (e.g., to impose certain types of priors).') ...
    }, ...
    'tv2d', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.'), ...
    arg({'useGPU','UseGPU'},false,[],'Use GPU acceleration. This is experimental and requires that UnLocBox is started with GPU support enabled.'), ...
    }, ...
    'tv3d', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.'), ...
    arg({'useGPU','UseGPU'},false,[],'Use GPU acceleration. This is experimental and requires that UnLocBox is started with GPU support enabled.'), ...
    }, ...
    'trace', { ...
    arg({'A','LinearOperator'},'@(x)x',[],'Linear transform. The norm applies to the linearly transformed feature vector. Either an expression that is evaluated in the workspace or a function handle.'), ...
    arg({'shape','Shape'},'',[],'Final feature shape. Allows to reshape the linearly transformed features into a matrix to apply matrix norms. If empty defaults to the shape of the original features.'), ...
    }}, 'Regularization param. Defines a param in the optimization problem; multiple types are supported and can be mixed freely.');

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
    arg({'regweights','TermWeights'},{{[]}},[],'Weights of the regularizers. This is a cell array of vectors of (relative) regularization parameters. Default is 1/N for N regularization terms. The cell array lists all possible assignments to search over.','type','expression','shape','row'), ...
    arg_sub({'solverOptions','SolverOptions'},{},{ ...
        arg({'maxit','MaxIterations'},1000,[],'Maximum number of iterations.'), ...
        arg({'rel_tol','RelativeTolerance'},2e-4,[],'Relative tolerance criterion. If the relative difference between two successive iterates is lower than this value the algorithm terminates.'),...
        arg({'abs_tol','AbsoluteTolerance'},0.000001,[],'Absolute tolerance criterion. If the objective function value falls below this the algorithm terminates.'), ...
        arg({'rho','CouplingParameter'},4,[],'Initial coupling parameter. For proximal algorithms this is the coupling strength between the terms between updates. Increasing this can improve the convergence speed but too strong values can prevent any convergence.'), ...
        arg({'rho_update','RhoUpdate'},true,[],'Update Rho. Whether to update rho dynamically according to 3.4.1 in [1]. Note, this can sometimes cause r_norm, s_norm to "blow up".'), ...
        arg({'rho_cutoff','RhoUpdateThreshold'},10.0,[],'Rho update threshold.','guru',true), ...
        arg({'rho_incr','RhoUpdateIncr'},2.0,[],'Rho update increment factor.','guru',true), ...
        arg({'rho_decr','RhoUpdateDecr'},2.0,[],'Rho update decrement factor.','guru',true), ...
        arg({'warmstart','Warmstart'},true,[],'Warm-start through regularization path. Enabling this is more efficient but convergence issues can be harder to trace down.') ...
    }, 'Controls the behavior of the ADMM optimization algorithm.'), ...
    arg_sub({'lbfgsOptions','LBFGSOptions'},{},{ ...
        arg({'MaxIter','MaxIterations'},10,[],'Maximum number of iterations.'), ...
        arg({'m','HessianHistory'},6,[],'LBFGS history length. The number of corrections to approximate the inverse hessian matrix.','guru',true), ...
        arg({'epsilon','Epsilon'},1e-4,[],'Tolerance criterion.  A minimization terminates when ||g|| < epsilon*max(1,||x||).'), ...
        arg({'past','DeltaDistance'},0,[],'Distance for delta-based convergence test.','guru',true), ...
        arg({'delta','Delta'},1e-5,[],'Delta for convergence test.','guru',true), ...
        arg({'linesearch','LineSearchAlgorithm'},'more_thuente',{'more_thuente','backtracking_armijo','backtracking_wolfe','backtracking_strong_wolfe'},'The line search algorithm.','guru',true), ...
        arg({'max_linesearch','MaxLineSearch'},40,[],'Maximum number of trials for the line search.','guru',true), ...
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
        arg({'lambdas','Lambdas'}, 2.^(3:-0.66:-8), [], 'Regulariation parameters. Controls the sparsity/simplicity of the result. Typically, this is an interval to scan, such as 2.^(10:-1:-15).'), ...
        arg({'nfolds','NumFolds'},5,[],'Cross-validation folds. The cross-validation is used to determine the best regularization parameter.'),...
        arg({'foldmargin','FoldMargin'},0,[],'Margin between folds. This is the number of trials omitted between training and test set.'), ...
        arg({'cvmetric','ParameterMetric'},'',{'','kld','nll','mcr','mae','mse','max','rms','bias','medse','auc','cond_entropy','cross_entropy','f_measure'},'Metric for Parameter Optimization. By default auto-determined; can be any of the ml_calcloss-supported metrics. In particular, auc is a good idea if the classification task is between highly imbalanced classes.') ...
        arg({'return_regpath','ReturnRegpath'}, true, [], 'Return the entire regularization path. This is for the best relative weighting of terms. If false, only the best model will be returned.'), ...
        arg({'return_reggrid','ReturnReggrid'}, false, [], 'Return the entire regularization grid. This also returns regularization paths for all other relative weightings.'), ...
        arg({'history_traces','HistoryTraces'}, false, [], 'Return history traces. If true, optimization history traces will be returned.'), ...
    }, 'Controls the search for the optimal regularization parameter.'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg_nogui({'shape','Shape'}, [], [], 'Reshaping for features. Allows to reshape (perhaps vectorized) features into a particular representation.'), ...
    arg({'verbosity','Verbosity'},1,[],'Diagnostic output level. Zero is off, 1 only shows cross-validation diagnostics, 2 shows solver diagnostics, 3 shows iteration diagnostics.'), ...
    arg({'continuous_targets','ContinuousTargets','Regression'}, false, [], 'Whether to use continuous targets. This allows to implement some kind of damped regression approach.'),...
    arg({'includebias','IncludeBias','bias'},true,[],'Include bias param. Also learns an unregularized bias param (strongly recommended for typical classification problems).'));

% find all target classes (if classification)
if iscell(targets)
    classes = unique([targets{:}]);
else
    classes = unique(targets);
end

if length(classes) > 2 && strcmp(loss,'logistic') && ~continuous_targets
    % in the multi-class case we use the voter for now (TODO: use softmax loss instead)
    model = ml_trainvote(trials, targets, '1v1', @ml_trainproximal, @ml_predictproximal, varargin{:});
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    % optionally allow nfolds to be a function of #trials
    if lambdaSearch.nfolds < 1
        lambdaSearch.nfolds = round(lambdaSearch.nfolds*length(targets)); end
    % sanitize some more inputs
    solverOptions.verbose = max(0,verbosity-1);    
    if isnumeric(regweights)
        regweights = {regweights}; end
    
    % lambdas need to be sorted in descending order for the warm-starting to work
    lambdaSearch.lambdas = sort(lambdaSearch.lambdas,'ascend');
    if strcmp(lambdaSearch.cvmetric,'mcr')
        lambdaSearch.cvmetric = ''; end
    
    % vectorize data if necessary & sanity-check explicit shape parameter if specified
    if ndims(trials) > 2 %#ok<ISMAT>
        featureshape = size(trials); featureshape = featureshape(1:end-1);
        if ~isempty(shape) && ~isequal(shape,featureshape)
            if prod(featureshape) == prod(shape)
                warning('You are specifying a shape property but also multidimensional features of a different shape; using the explicit shape parameter.');
                featureshape = shape;
            else
                error('You are specifying a shape property but also features with an incompatible number of elements. Please correct.');
            end
        end
        trials = double(reshape(trials,[],size(trials,ndims(trials)))');
        vectorize_trials = true;
    else
        if ~isempty(shape)
            if size(shape,1) > 1
                if all(all(bsxfun(@eq,shape(1,:),shape)))
                    featureshape = [shape(1,1),shape(1,2),size(shape,1)];
                    if prod(featureshape) == size(trials,2)
                        % the reason is that it is much more efficient to operate on a dense 3d array than a very sparse 2d array
                        warn_once('This method will by convention reshape block-diagonalized feature matrices with identical blocks into a 3d tensor. This warning will only come up once.');
                    else
                        error('Your shape parameter has a different number of features than your data.');
                    end
                else
                    % we don't implement block-diagonalization in here
                    error('This method does not handle implicitly block-diagonal features; please either reformulate in tensor form or pass a large sparse data matrix (pre-blockdiagonalized). Note that the tensor form is likely several times faster.');
                end
            elseif prod(shape) ~= size(trials,2)
                error('Your shape parameter has a different number of features than data.');
            else
                featureshape = shape;
            end
        else
            featureshape = [size(trials,2),1];
        end
        vectorize_trials = false;
    end
    
    % optionally scale the data
    sc_info = hlp_findscaling(trials,scaling);
    trials = hlp_applyscaling(trials,sc_info);
    
    % optionally remap target labels to -1,+1
    if strcmp(loss,'logistic') && length(classes) == 2 && ~continuous_targets
        targets(targets==classes(1)) = -1;
        targets(targets==classes(2)) = +1;
    end
        
    % learn a sequence of models across the given lambda's, on all the data (i.e. the regularization path)
    if verbosity
        disp('Running optimization...'); end
    
    % run a cross-validation to score the lambdas and regweights
    loss_mean = cell(1,length(regweights));
    predictions = repmat({zeros(length(targets),length(lambdaSearch.lambdas))},1,length(regweights)); % predictions{r}(t,l) is the prediction for regweight combination #r, trial #t and lambda #l
    reptargets = repmat(targets,1,length(lambdaSearch.lambdas));
    foldid = 1+floor((0:length(targets)-1)/length(targets)*lambdaSearch.nfolds);
    % for each fold...
    for i = 1:lambdaSearch.nfolds
        if verbosity
            disp(['Fitting fold # ' num2str(i) ' of ' num2str(lambdaSearch.nfolds)]); end

        % determine training and test sets
        which = foldid==i;
        trainids = ~which;
        whichpos = find(which);
        for j=1:lambdaSearch.foldmargin
            trainids(max(1,whichpos-j)) = false;
            trainids(min(length(which),whichpos+j)) = false;
        end
        testset = [trials(which,:) ones(length(whichpos),double(includebias))];

        % for each relative regularization term weighting...
        for w = 1:length(regweights)

            % get regularization path
            [model_seq{w},history_seq{w}] = hlp_diskcache('predictivemodels',@solve_regularization_path,trials(trainids,:),targets(trainids),lambdaSearch.lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,regweights{w},featureshape);  %#ok<NASGU>
            
            % calc test-set predictions for each model
            for m=1:length(model_seq{w})
                predictions{w}(which,m) = (testset*model_seq{w}{m}(:))'; end
            if strcmp(loss,'logistic')
                predictions{w}(which,:) = 2*(1 ./ (1 + exp(-predictions{w}(which,:))))-1; end
            
            % evaluate losses on the fold
            if isempty(lambdaSearch.cvmetric)
                if strcmp(loss,'logistic')
                    loss_mean{w}(i,:) = mean(reptargets(which,:) ~= sign(predictions{w}(which,:)));
                else
                    loss_mean{w}(i,:) = mean((reptargets(which,:) - predictions{w}(which,:)).^2);
                end
            else
                for r=1:length(lambdaSearch.lambdas)
                    loss_mean{w}(i,r) = ml_calcloss(lambdaSearch.cvmetric,reptargets(which,r),predictions{w}(which,r)); end
            end
        end
        1;
    end

    % find best lambdas & regweight combination
    for k=1:length(regweights)
        % average over folds
        loss_mean{k} = mean(loss_mean{k},1);
        % if there are several minima, choose largest lambda of the smallest cvm
        best_lambda{k} = max(lambdaSearch.lambdas(loss_mean{k} <= min(loss_mean{k})));
        best_loss{k} = min(loss_mean{k});
    end
    % pick regweights and lambda at best loss
    [dummy,best_regidx] = min([best_loss{:}]); %#ok<ASGLU>
    best_regweights = regweights{best_regidx};
    lambda_min = best_lambda{best_regidx};

    % pick the model at the minimum...
    if lambdaSearch.return_regpath
        [regpath,history] = hlp_diskcache('predictivemodels',@solve_regularization_path,trials,targets,lambdaSearch.lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,best_regweights,featureshape);
        model.regularization_path = regpath;
        model.regularization_loss = loss_mean{best_regidx};
        model.w = regpath{find(lambdaSearch.lambdas == lambda_min,1)};
    else
        [tmp,history] = hlp_diskcache('predictivemodels',@solve_regularization_path,trials,targets,lambdaSearch.lambdas(find(lambdaSearch.lambdas==lambda_min,1)),loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizers,best_regweights,featureshape);
        model.w = tmp{1};
    end
    model.balance_losses = loss_mean;
    if lambdaSearch.history_traces
        model.fold_history = history_seq;
        model.regularization_history = history;
    end
    if lambdaSearch.return_reggrid
        model.regularization_grid = model_seq; end
    model.classes = classes;
    model.continuous_targets = continuous_targets;
    model.includebias = includebias;
    model.vectorize_trials = vectorize_trials;
    model.featureshape = featureshape;
    model.sc_info = sc_info;
    model.loss = loss;
end



% learn the regularization path
function [regpath,fullhist] = solve_regularization_path(A,y,lambdas,loss,includebias,verbosity,solverOptions,lbfgsOptions,regularizersArg,regweights,featureshape)
% solve_regularization_path_version<1.0.0>
fullhist = {};

[m,n] = size(A); % m trials, n features
w = zeros(n+double(includebias),1);

% derive the design matrix A & label vector y from the trials...
A = double(A);
y = double(y(:));

% set up the data-dependent loss function to use
switch loss
    case 'logistic'
        C = [bsxfun(@times,-y,A) -y];
        if lbfgsOptions.useGPU
            try
                C = gpuArray(C); 
            catch e
                disp_once(['Could not enable GPU support: ' e.message]);
            end
        end
        Cp = C';
        lossfunc.prox = @(x,gamma,x0) prox_logistic(C,Cp,x,gamma,x0,m,hlp_struct2varargin(lbfgsOptions));
        lossfunc.eval = @(x,lambda) lambda*obj_logistic(C,x,n);
    case 'squared'
        % append a bias to the design matrix
        Ao = [A ones(size(A,1),1)];
        [mm,nn] = size(Ao);
        % choose the right prox operator
        if solverOptions.rho_update
            lossfunc.prox = @(x,gamma,x0) prox_squared_iter(Ao,y,1/gamma,x,zeros(size(x)),nn,x0);
        else            
            Atb = Ao'*y;
            [L,U] = factor(Ao,solverOptions.rho);
            lossfunc.prox = @(x,gamma,x0) prox_squared_factored(Ao,Atb,L,U,solverOptions.rho,x,zeros(size(x)),mm,nn);
        end
        lossfunc.eval = @(x,lambda) lambda*obj_squared(Ao,y,x);
    case 'hyperbolic-secant'
        % lossfunc.prox = @(x,gamma,x0) prox_hs(y,x,gamma,x0,hlp_struct2varargin(lbfgsOptions));
        % lossfunc.eval = @(x,lambda) lambda*obj_hs(x,b);
        error('Hyperbolic-secant loss is not yet implemented.');
    otherwise
        error('Unsupported loss function.');
end
lossfunc.y0 = [];
lossfunc = @(lambda)setfield(setfield(lossfunc,'prox',@(x,gamma,x0)lossfunc.prox(x,gamma*lambda,x0)),'eval',@(x)lossfunc.eval(x,lambda)); %#ok<SFLD>

% ensure that regularizers is a cell array
regularizers = {};
if isstruct(regularizersArg)
    for k=1:length(fieldnames(regularizersArg))
        if isfield(regularizersArg,['term' num2str(k)])
            regularizers{end+1} = regularizersArg.(['term' num2str(k)]); end
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
                [a,b,c,d,e,f,g,h] = size(reshape(w(1:n),featureshape)); %#ok<ASGLU>
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
            rA(reshape(w(1:n),featureshape));
            param.A = @(x)rA(reshape(x(1:n),featureshape));
        catch
            try
                param.A(w(1:n));
                param.A = @(x)rA(x(1:n));
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
            [a,b,c,d,e,f,g,h] = size(reshape(w(1:n),featureshape)); %#ok<ASGLU>
            param.shape = eval(param.shape);
        end
        
        % set the A matrix for future reference
        if isfield(param,'A')
            rA = param.A;
        else
            rA = @(x)x(1:n);
        end
        
        % move the A parameter into regfunc.L (handled by ADMM)
        if isfield(param,'A') && ~(strcmp(param.arg_selection,'l2') && param.nonorthogonal_transform)
            regfunc.L = param.A;
            % remove fields from param
            param = rmfield(param,'A');
            if isfield(param,'At')
                param = rmfield(param,'At'); end
        else
            regfunc.L = @(x)x(1:n);
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
                    A = operator_to_matrix(param.A,n);
                    regfunc.prox = @(x,gamma,x0) prox_squared_iter(A,param.y,1/gamma,x,zeros(size(x)),n,x0);
                    regfunc.eval = @(x,lambda) lambda*obj_squared(A,param.y,x(1:n));
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
        regfuncs{end+1} = @(lambda) setfield(setfield(regfunc,'prox',@(x,gamma,x0)regfunc.prox(x,gamma*lambda,x0)),'eval',@(x)regfunc.eval(x,lambda)); %#ok<SFLD>
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
regpath = cell(1,length(lambdas));
for k =1:length(lambdas)
    
    % set up parameters
    weights = [1,lambdas(k)*regweights];
    lossfunc = lossfunc(1);
    lossfunc.L = [];
    if ~isempty(y0)
        lossfunc.y0 = y0{1}; end
    lossfunc.param = struct();
    for r = 1:length(regfuncs)
        tmpregfuncs(r) = regfuncs{r}(weights(1+r));
        if ~isempty(y0)
            tmpregfuncs(r).y0 = y0{1+r}; end
    end
    
    % solve
    t0 = tic;
    if verbosity
        fprintf('  scanning lambda = %f (%i/%i)...',lambdas(k),k,length(lambdas)); end
    if solverOptions.warmstart
        [w,y0,rho,hist] = hlp_diskcache('temporary',@consensus_admm,w,[lossfunc tmpregfuncs],solverOptions); %#ok<ASGLU>
    else
        [w,y0dummy,rho,hist] = hlp_diskcache('temporary',@consensus_admm,zeros(size(w)),[lossfunc tmpregfuncs],solverOptions); %#ok<ASGLU>
    end
    if verbosity
        fprintf(' %i iters; t = %.1fs\n',length(hist.objval),toc(t0)); end
    
    % store
    regpath{k} = w;
    fullhist{k} = hist;
end

% --- logistic loss code ---

function [val,grad] = obj_proxlogistic(x,C,Cp,z,gamma,m)
% objective function for the logistic loss proximity operator (effectively l2-regularized logreg)
ecx = exp(C*x);
val = (1/2)*sum((x-z).^2) + gamma/m*gather(sum(log1p(ecx)));
if ~isfinite(val)
    ecx(~isfinite(ecx(:))) = 2.^50;
    val = (1/2)*sum((x-z).^2) + gamma/m*gather(sum(log1p(ecx)));
end
grad = (x - z) + (gamma/m)*gather(Cp*(ecx./(1+ecx)));

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


function M = operator_to_matrix(op,n)
persistent operator_cache;
opid = ['x' hlp_cryptohash({op,n})];
try
    M = operator_cache.(opid);
catch
    % convert a linear operator function to a matrix (given the input dimensionality)
    % this works by going through the canonical basis vectors of the space and projecting them
    % one-by-one through the operator
    fprintf('Evaluating and caching operator...');
    M = hlp_diskcache('general',@operator_to_matrix_cached,op,n);
    operator_cache.(opid) = M;
    fprintf('done.\n');
end

function M = operator_to_matrix_cached(op,n)
vec = @(x)x(:);
M = sparse([]);
w = zeros(n,1);
for c=[n 1:n-1]
    v = w; v(c) = 1;
    M(:,c) = vec(op(v));
end
