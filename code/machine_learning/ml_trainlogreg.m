function model = ml_trainlogreg(varargin)
% Learn a linear probabilistic predictive model by logistic regression.
% Model = ml_trainlogreg(Trials, Targets, Lambda, Options...)
%
% Logistic regression is the simplest fully probabilistic classifier in the toolbox, and can be
% understood as a linear projection followed by a mapping through a sigmoid (0-1) function, which
% allows to translate the input space smoothly into probabilities. Logreg is available in several
% variants, first the three variational Bayes (vb) ones, which are not dependent on a regularization
% parameter [1]. 'vb' and 'vb-iter' are about as versatile as LDA, except that they are fully
% probabilistic and not as outlier-prone; they differ mostly in numerical robustness. 'vb-ard' uses
% automatic relevance determination [2] to obtain a mapping that is sparse in the features, i.e.
% only a subset of the features is used for prediction. The Bayesian implementations compute not
% only a point estimate of the weights, but a full posterior distribution (gaussian), which enables
% more advanced post-processing of the model and outputs.
%
% The l1/l2 variants are computed via convex optimization. They might generally be somewhat more
% robust in extreme corner cases than the variational variants. Both require that a regularization
% parameter is specified, which can in simple cases be left at the default, but which is ideally
% optimized using parameter search, especially in the presence of many features and/or
% hard-to-handle distributions (more important for l1 than l2). This is computationally expensive,
% however. The 'l2' variant is the general-purpose one and in practice very similar to the
% vb/vb-iter variants. The 'l1' variant uses an l1-norm on the weights as regularizer and therefore
% gives a model that is sparse in the features (like vb-ard).
%
% A particularly fast implementation that is equivalent to l1 is LARS (least-angle regression) [7], 
% which is generally to be preferred for sparse results, when speed is an issue and the binary is 
% executable.
%
% Logistic regression is a very powerful, robust and versatile classifier for biosignals, and also
% very fast (when not regularized) - it is worth trying it on almost every problem that is hoped to
% be (more or less) linearly separable. Sparse logistic regression is a special case which allows to
% select features (out of exponentially many irrelevant ones) [3,4] on top of doing classification.
% EEG-derived feature are however in most cases not sparse, so that l1 or vb-ard may give worse
% results. Exceptions are features derived from independent components and other sparsity-inducing
% feature extractors. Some signal bases (e.g. wavelets) may be half-way in between sparse and
% non-sparse.
%
% The Bayesian variants are implemented using the bayes_logit toolbox by Jan Drugowitsch, and the
% regularized variants are implemented using the LIBLINEAR package [5] or CVX [6] as a fall-back.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%                  may also contain weights, in the form {Targets,Weights}; supported by l1 & l2 variants
%
%   Lambda       : regularization parameter, for l1/l2 logistic regression variant (default: 1)
%                  a comprehensive search interval would be 2.^(-6:0.5:10)
%
%   Options  : optional name-value parameters to control the training details:
%              'variant': one of several logreg variants:
%                        'vb': variational Bayes (Jakkola&Jordan, 2000), using a joint uninformative
%                              prior on the weights, but possibly numerically unstable in corner cases (default)
%                        'vb-iter': variational Bayes (Jakkola&Jordan, 2000), iterative, using a per-weight prior,
%                                   more or less uninformative, but numerically more robust than the 'vb' variant
%                        'vb-ard': variational Bayes (Jakkola&Jordan, 2000), implements sparse logistic regression
%				                   (using automated relevance determination)
%                        'lars' : using least-angle regression (efficient version of 'l1', using GLMNET)
%                                 note: this can also be specified as a cell array of the form 
%                                 {'lars', name,value, name,value, ...} where the names and values 
%                                 specify custom options for the glmnet solver.
%                        'l1' : sparse logistic regression using l1 regularization, using either 
%                               LIBLINEAR or CVX
%                        'l2' : logistic regression using l2 regularization using LIBLINEAR or CVX
%                       
%
%              'eps'    : desired accuracy (default: [] - the default of the respective library, currently only supported for l1/l2)
%
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std')
%
%              'regression' : do binomially-distributed regression as opposed to binary classification (default: false)
%
% Out:
%   Model   : a linear model;
%             w is the linear weights, b is the bias;
%             V is the covariance matrix of a posterior distribution around the weights (jointly over w and b)
%               (only available for the three 'vb' modes)
%             classes indicates the class labels which the model predicts
%             additional parameters determine a posterior distribution over the weights
%
% Examples:
%   % learn a logistic regression model using a fast Bayesian method
%   model = ml_trainlogreg(trials,targets)
%
%   % as before, but using a numerically more robust (yet slower) approach
%   model = ml_trainlogreg(trials,targets,[],'variant','vb-iter')
%
%   % learn a sparse Bayesian logistic regression model 
%   model = ml_trainlogreg(trials,targets,[],'variant','vb-ard')
%
%   % learn a sparse regularized logistic regression model using LARS
%   model = ml_trainlogreg(trials,targets,[],'variant','lars')
%
%   % as before, but use a different number of inner cross-validation folds for LARS.
%   model = ml_trainlogreg(trials,targets,[],'variant',{'lars','NumFolds',10})
%
%   % learn a sparse regularized logistic regression model using a relatively less efficient approach
%   model = utl_searchmodel({trials,targets},'args',{{'logreg',search(2.^(-6:1:10)),'variant','l1'}})
%
%   % learn a regularized logistic regression model (using cross-validation)
%   model = utl_searchmodel({trials,targets},'args',{{'logreg',search(2.^(-6:1:10)),'variant','l1'}})
%
% See also:
%   ml_predictlogreg
%
% References:
%   [1] Jaakkola, T. S., and Jordan, M. I. "A variational approach to bayesian logistic regression models and their extensions."
%       In Proceedings of the Sixth International Workshop on Artificial Intelligence and Statistics (1997).
%   [2] Wipf, D., Nagarajan, S., Platt, J., Koller, D., Singer, Y., and Roweis, S. "A New View of Automatic Relevance Determination."
%       MIT Press, 2008, pp. 1632, 1625.
%   [3] P. Zhao and B. Yu. "On model selection consistency of Lasso."
%       JMLR, 7:2541–2563, 2006.
%   [4] M. J. Wainwright. "Sharp thresholds for noisy and high-dimensional recovery of sparsity using l1-constrained quadratic programming."
%      Technical Report 709, Dpt. of Statistics, UC Berkeley, 2006.
%   [5] R.-E. Fan, K.-W. Chang, C.-J. Hsieh, X.-R. Wang, and C.-J. Lin. "LIBLINEAR: A library for large linear classification"
%       Journal of Machine Learning Research 9(2008), 1871-1874.
%   [6] M. Grant and S. Boyd. "Graph implementations for nonsmooth convex programs"
%       Recent Advances in Learning and Control (a tribute to M. Vidyasagar), V. Blondel, S. Boyd, and H. Kimura, editors, pages 95-110,
%       Lecture Notes in Control and Information Sciences, Springer, 2008.
%   [7] Efron B., Hastie T., Johnstone I., and Tibshirani R., "Least Angle Regression",
%       Annals of Statistics 32(2), 407-499. 2004.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'lam','Lambda','lambda'},  1, [0 2^-7 2^15 Inf], 'Regularization parameter for l1/l2. A comprehensive search interval would be 2.^(-6:0.5:10).'), ...
    arg_subswitch({'variant','Variant'},'vb', ...
        {'vb',{},'vb-iter',{},'vb-ard',{},'l1',{},'l2',{}, ...
         'lars',{arg({'nfolds','NumFolds'},5,[0 Inf],'Cross-validation folds. The cross-validation is used to determine the best regularization parameter. If this is smaller than 1, it is taken as a fraction of the number of trials.'),...
                  arg({'foldmargin','FoldMargin'},0,[0 0 10 Inf],'Margin between folds. This is the number of trials omitted between training and test set.'), ...
                  arg({'alpha','ElasticMixing'},1,[0.01 1],'ElasticNet mixing parameter. The default is the lasso penalty.'), ...
                  arg({'nlambda','NumLambdas'},100,uint32([1 10 500 5000]),'Number of lambdas. Number of lambda (regularization) parameter values to consider.','guru',true), ...
                  arg({'penalty_factor','PenaltyFactor'},[],[],'Penalty per feature. Allows for selective shrinkage and the like.'), ...
                  arg({'lambda_min','MinLambda'},0,[0 Inf],'Minimum lambda. Smallest value for lambda, as a fraction of lambda_max, the (data derived) entry value (i.e., the smallest value for which all coefficients are zero). By default 0.0001, or 0.05 if #observations < #trials.','guru',true), ...
                  arg({'thresh','ConvergenceThreshold'},[],[0 Inf],'Convergence threshold. Each inner coordinate-descent loop continues until the relative change in any coefficient is less than this (if specified).','guru',true), ...
                  arg({'dfmax','MaxVariables'},0,uint32([0 1000000000]),'Max #of variables. If specified, this constraints the maximum number of variables in the model.','guru',true), ...
                  arg({'pmax','MaxNonzeroes'},0,uint32([0 1000000000]),'Max #of non-zeroes. If specified, this constraints the maximum number of non-zero variables in the model.','guru',true), ...
                  arg({'maxit','MaxIterations'},300,uint32([1 50 1000 10000]),'Max iterations. The maximum number of iterations, if the convergence threshold is not reached.'), ...
                  arg({'HessianExtract','ExtractHessian'},false,[],'Extract Hessian. If false (the default), an upper-bound approximation is made to the hessian, which is not recalculated at each outer loop..','guru',true), ...
                  arg({'verbose','VerboseOutput'},false,[],'Verbose outputs. Whether to display verbose console outputs.'), ...
                  arg({'customloss','CustomLoss'},'default',{'default','auto','kld','nll','mcr','mae','mse','max','rms','bias','medse','auc','cond_entropy','cross_entropy','f_measure'},'Search lambda according to a custom loss. Note that some loss types may not be applicable here and give errors.')} ...
        }, 'Variant to use. Variational Bayes methods: using a joint uniformative prior (possibly unstable), or using a per-weight prior (iter), or using a sparse prior (Automatic Relevance Determination). Regularized methods: lars for fast sparse logistic regression, l1 for slow sparse logistic regression, or l2 for not necessarily sparse results.'),...
    arg({'epsi','Epsilon','eps'}, [], [], 'Desired residual error. Currently only supported for the l1/l2 variants.'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'),...
    arg({'continuous_targets','ContinuousTargets','Regression'}, false, [], 'Whether to use continuous targets. This allows to implement some kind of damped regression approach.'),...
    arg({'votingScheme','VotingScheme'},'1v1',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'), ...
    arg({'fallback','UseFallback'}, false, [], 'Use CVX fallback. The CVX fallback does not depend on binaries and should yield consistent results across all platforms.'));

if is_search(lam)
    lam = 1; end
if ~isnumeric(lam)
    error(['Lambdas is expected to be numeric, but was assigned: ' hlp_tostring(lam)]); end

% obtain weights
if iscell(targets) %#ok<*NODEF>
    [targets,weights] = deal(targets{:});
    weights = weights/sum(weights);
else
    weights = [];
end

% we have to use CVX for logreg, if the appropriate LIBLINEAR version is missing
need_fallback = fallback || isempty(which('llwtrain'));
use_cvx = any(strcmp(variant.arg_selection,{'l1','l2'})) && need_fallback;

% identify and remap the classes
classes = unique(targets);
if length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.'); end

% optionally remap target values to the proper class indices
if ~continuous_targets
    if ~isempty(strfind(variant.arg_selection,'vb')) || use_cvx
        % these are two-class classifiers
        if length(classes) > 2
            % so we need to vote in this case
            model = ml_trainvote(trials,targets,votingScheme,@ml_trainlogreg,@ml_predictlogreg,varargin{:});
            return;
        end
        % remap target labels to -1,+1
        targets(targets==classes(1)) = -1;
        targets(targets==classes(2)) = +1;
    elseif strcmp(variant.arg_selection,'lars')
        % remap target labels to 1..k
        targ = targets;
        for c=1:length(classes)
            targets(targ==classes(c)) = c; end
    else
        % remap target labels to 0..k-1
        targ = targets;
        for c=1:length(classes)
            targets(targ==classes(c)) = c-1; end
    end
end


% scale the data
sc_info = hlp_findscaling(trials,scaling);
trials = hlp_applyscaling(trials,sc_info);

% add bias variable to training data for the VB's
if ~isempty(strfind(variant.arg_selection,'vb'))
    warning off Bayes:maxIter
    trials = [trials ones(size(trials,1),1)];
end

switch variant.arg_selection
    case 'vb'
        [model.w, model.V, model.invV, model.logdetV] = hlp_diskcache('predictivemodels',@bayes_logit_fit,trials,targets);
    case 'vb-iter'
        [model.w, model.V, model.invV, model.logdetV] = hlp_diskcache('predictivemodels',@bayes_logit_fit_iter,trials,targets);
    case 'vb-ard'
        [model.w, model.V, model.invV, model.logdetV] = hlp_diskcache('predictivemodels',@bayes_logit_fit_ard,trials,targets);
    case 'lars'
        % build options structure
        if isempty(variant.thresh)
            variant = rmfield(variant,'thresh'); end
        if isfield(variant,'arg_direct')
            variant = rmfield(variant,'arg_direct'); end
        glmopts = variant;
        % add weights, if given
        if ~isempty(weights)
            glmopts.weights = weights; end
        if variant.nfolds < 1
            variant.nfolds = round(size(trials,1)*variant.nfolds); end
        if ~continuous_targets
            family = 'multinomial';
        else
            family = 'gaussian'; % note: glmnet can't handle logistic regression with continuous targets, so need to revert to Gaussian
        end
        % run
        foldid = 1+floor((0:length(targets)-1)/length(targets)*variant.nfolds);
        model.CVerr = cvglmnetMulticlass(double(trials),double(targets),variant.nfolds,variant.foldmargin,foldid,'response',family,glmnetSet(glmopts),variant.verbose,variant.customloss);
        try
            % assign reasonable weights for visualization
            model.w = model.CVerr.glmnet_object.beta{2}(:,model.CVerr.glmnet_object.lambda == model.CVerr.lambda_min);
        catch
            disp_once('ml_trainlogreg: Could not assign weights for visualization.');
        end
    case 'l1'
        lam = 1/lam; % lambda is parameterized differently in LIBLINEAR, we adapt to that
        try
            if use_cvx
                error('fall back'); end
            if isempty(weights)
                weights = ones(size(targets)); end
            % LIBLINEAR-weights 1.6+
            model = hlp_diskcache('predictivemodels',@llwtrain,double(weights),double(targets),sparse(double(trials)),[...
                '-s 6 -c ' num2str(lam) ' -B 1' fastif(~need_fallback,' -q ','') ...
                fastif(~isempty(epsi),[' -e ' num2str(epsi)],'')]);
        catch
            % fallback
            if ~isempty(weights)
                cvx_begin
                    variables W(size(trials,2)) b
                    minimize(sum(weights.*log(1+exp(-targets.*(trials*W+b))))*lam + norm(W,1))
                cvx_end
            else
                cvx_begin
                    variables W(size(trials,2)) b
                    minimize(sum(log(1+exp(-targets .* (trials*W+b))))*lam + norm(W,1))
                cvx_end                
            end
            model = struct('W',W,'b',b);
        end
    case 'l2'
        lam = 1/lam;
        try
            if use_cvx
                error('fall back'); end
            if isempty(weights)
                weights = ones(size(targets)); end
            % LIBLINEAR-weights 1.6+
            model = hlp_diskcache('predictivemodels',@llwtrain,double(weights),double(targets),sparse(double(trials)),[...
                '-s 0 -c ' num2str(lam) ' -B 1' fastif(~need_fallback,' -q ','') ...
                fastif(~isempty(epsi),[' -e ' num2str(epsi)],'')]);
        catch
            if ~isempty(weights)
                % fallback
                cvx_begin
                    variables W(size(trials,2)) b
                    minimize(sum(weights.*log(1+exp(-targets.*(trials*W+b))))*lam + norm(W,2))
                cvx_end
            else
                cvx_begin
                    variables W(size(trials,2)) b
                    minimize(sum(log(1+exp(-targets .* (trials*W+b))))*lam + norm(W,2))
                cvx_end
            end
            model = struct('W',W,'b',b);
        end
end

model.classes = classes;
model.sc_info = sc_info;
model.variant = variant.arg_selection;
model.use_cvx = use_cvx;
model.continuous_targets = continuous_targets;