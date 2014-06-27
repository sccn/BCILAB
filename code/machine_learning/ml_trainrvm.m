function model = ml_trainrvm(varargin)
% Learn a probabilistic (non-)linear model, via the Relevance Vector Machine.
% Model = ml_trainrvm(Trials, Targets, Options...)
%
% The Relevance Vector Machine [1] is a Bayesian equivalent to the popular Support Vector Machines
% (SVMs) [2]. RVMs can be used for general-purpose linear or non-linear (kernelized) classification
% or regression, and produce state-of-the-art results in most cases. Various kernels are supplied,
% where the rbf kernel is usually the best initial choice. In the non-linear case, the kernel
% scaling parameter should be found via parameter search, which can be time-consuming. In contrast
% to SVMs, Relevance Vector Machines give probablistic outputs, which can be practical when multiple
% uncertain predictions are to be fused, etc. In the future, an implementation of RVMs using convex
% optimization will be provided, which is assumed to give better optimality guarantees than the
% current implementation [3,4,5].
%
% RVMs (as well as kernel SVMs) are the most versatile general-purpose classifiers currently
% available in the toolbox, and are a good default choice if very little is known about the data.
% For best results, care must be taken to search over the respective regularization parameter(s)
% appropriately. If more is known about the structure of the data (e.g. that it should be linearly
% separable, or that that sparsity can be exploited across features), specialized methods may be
% more appropriate (and give similar results faster or reach better performance by overfitting less
% strongly).
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Options  : optional name-value parameters to control the training details:
%              'ptype': problem type: 'classification' (default) or 'regression'
%
%              'kernel': one of several kernel types:
%                         * 'linear':   Linear
%                         * 'rbf':      Gaussian / radial basis functions (default)
%                         * 'laplace':  Laplacian
%                         * 'poly':		Polynomial
%                         * 'cauchy':	Cauchy
%
%              'gamma': scaling parameter of the kernel (for regularization); if multiple values
%                       are given, the optimal gamma will be searched via evidence maximization
%                       default: 2.^(-16:0.2:10)
%
%              'degree': degree of the polynomial kernel, if used (default: 3)
%
%              'bias': whether to add a bias to the data (default: 1)
%
%              misc options:
%              'iterations': Number of interations to run for
%              'time':       time limit to run for, e.g. '1.5 hours', '30 minutes', '1 second'
%              'fixednoise': whether the gaussian noise is to be fixed 0/1
%              'beta':       (Gaussian) noise precision (inverse variance)
%              'noisestd':   (Gaussian) noise standard deviation
%              'scaling':    pre-scaling of the data (see hlp_findscaling for options) (default: 'std') 
%              'diagnosticlevel':   verbosity level, 0-4
%
% Out:
%   Model   : the computed model...
%             classes indicates the class labels which the model predicts
%             additional parameters determine a posterior distribution over the weights
%
% Examples:
%   % learn a standard Relevance Vector Machine classifier
%   model = ml_trainrvm(trials,targets)
%
%   % as before, but this time use a regression approach
%   model = ml_trainrvm(trials,targets,'ptype','regression')
%
%   % use a Laplacian kernel 
%   model = ml_trainrvm(trials,targets,'kernel','laplace')
%
%   % find the optimal kernel scale using parameter search
%   model = utl_searchmodel({trials,targets},'args',{{'rvm','gamma',seach(2.^(-16:2:4)))
%
%   
% See also:
%   ml_predictrvm, SparseBayes
%
% References:
%  [1] Vladimir Vapnik. "The Nature of Statistical Learning Theory." 
%      Springer-Verlag, 1995
%  [2] Michael E. Tipping and Alex Smola, "Sparse Bayesian Learning and the Relevance Vector Machine". 
%      Journal of Machine Learning Research 1: 211?244. (2001)
%  [3] Michael E. Tipping and A. C. Faul. "Fast marginal likelihood maximisation for sparse Bayesian models."
%      In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence and Statistics, Key West, FL, Jan 3-6 (2003)
%  [4] David P. Wipf and Srikantan Nagarajan, "A New View of Automatic Relevance Determination,"
%      In J.C. Platt, D. Koller, Y. Singer, and S. Roweis, editors, Advances in Neural Information Processing Systems 20, MIT Press, 2008.
%  [5] David P. Wipf and Srikantan Nagarajan, "Sparse Estimation Using General Likelihoods and Non-Factorial Priors," 
%      In Advances in Neural Information Processing Systems 22, 2009.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-06

opts = arg_define([0 2],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'ptype','Type'}, 'classification', {'classification','regression'}, 'Type of problem to solve.','cat','Core Parameters'), ...
    arg({'kernel','Kernel'}, 'rbf', {'linear','rbf','laplace','poly','cauchy'}, 'Kernel type. Linear, or Non-linear kernel types: Radial Basis Functions (general-purpose), Laplace (sparse), Polynomial (rarely preferred), and Cauchy (slightly experimental).','cat','Core Parameters'), ...
    arg({'gammap','KernelScale','gamma'}, 2.^(-16:0.5:10), [0 2^-20 2^10 Inf], 'Scaling of the kernel functions. Should match the size of structures in the data. A reasonable range is 2.^(-16:2:4).','cat','Core Parameters','shape','row'), ...
    arg({'polydegree','PolyDegree','degree'}, 3, uint32([1 100]), 'Degree of the polynomial kernel, if chosen.','cat','Core Parameters'), ...
    arg({'bias','Bias'}, true, [], 'Include a bias term in the model.','cat','Core Parameters'), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.','cat','Core Parameters'), ...
    ...
    arg({'iterations','MaxIterations'}, 100, uint32([1 100000]), 'Number of iterations to run.','cat','Miscellaneous'), ...
    arg({'time','MaxTime'}, '1000 seconds', [], 'Maximum time to run. Can use ''seconds'', ''minutes'', ''hours'' in the string.','cat','Miscellaneous'), ...
    arg({'fixednoise','NoiseFixed'}, false, [], 'Keep the Gaussian noise estimate fixed.','cat','Miscellaneous'), ...
    arg({'noiseinvvar','NoiseInvVariance','beta'}, [], [], 'Inverse variance of the Gaussian noise term.','cat','Miscellaneous','shape','scalar'), ...
    arg({'noisestd','NoiseVariance'}, [], [], 'Variance of the Gaussian noise term.','cat','Miscellaneous','shape','scalar'), ...
    arg({'diagnosticlevel','Verbosity'}, 'none', {'none','minimal','low','medium','high','ultra'}, 'Verbosity level.','cat','Miscellaneous'),...
    arg({'votingScheme','VotingScheme'},'1vR',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'), ...    
    arg({'monitor','DisplayInterval'}, uint32(0), [], 'Iterations between diagnostic outputs.','cat','Miscellaneous'));

arg_toworkspace(opts);

if is_search(gammap)
    gammap = 0.3; end

% pre-process arguments
ptype = hlp_rewrite(ptype,'classification','c','regression','r'); %#ok<*NODEF>
likelihood = hlp_rewrite(ptype,'c','bernoulli','r','gaussian'); 
args1 = [hlp_struct2varargin(opts,'restrict',{'iterations','time','monitor','fixednoise','freebasis','callback','callbackdata'}),{'diagnosticlevel',hlp_rewrite(opts.diagnosticlevel,'minimal','none')}];
args2 = hlp_struct2varargin(opts,'restrict',{'beta','noisestd','relevant','weights','alpha'},'rewrite',{'beta','noiseinvvar'});

% remap targets for classification
if strcmp(ptype,'c')
    classes = unique(targets);
    if length(classes) > 2
        % multiclass case: use the voter
        model = ml_trainvote(trials,targets,votingScheme,@ml_trainrvm,@ml_predictrvm,varargin{:});
        return;
    elseif length(classes) == 1
        error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
    end
    % remap target labels to 0/1
    targets(targets==classes(1)) = 0;
    targets(targets==classes(2)) = 1;
else
    classes = [];
end

% prescale the data
sc_info = hlp_findscaling(trials,scaling);
trials = hlp_applyscaling(trials,sc_info);
basis = trials;

if length(gammap)>1
    if ~strcmp(opts.diagnosticlevel,'none')
        disp('Now optimizing gamma parameter using evidence maximization...'); end
    % optimize gamma parameter using marginal log-likelihood    
    bestgam = NaN;      % best gamma parameter so far
    likelihoods = nan(length(gammap),1);
    bestlike = -Inf;    % best likelihood so far
    for k=1:length(gammap)
        gam = gammap(k);
        if ~strcmp(opts.diagnosticlevel,'none')
            fprintf('gamma=%.3f\n',gam); end
        % kernelize the data and add bias
        ktrials = utl_kernelize(trials,basis,kernel,gam,polydegree);
        ktrials = quickif(bias,[ones(size(ktrials,1),1) ktrials],ktrials);
        % run the RVM
        [param, hyperparam, diag] = hlp_diskcache('predictivemodels',@SparseBayes,likelihood,ktrials,targets,SB2_UserOptions(args1{:}),SB2_ParameterSettings(args2{:})); %#ok<ASGLU>        
        if ~isempty(diag.Likelihood)
            likelihoods(k) = diag.Likelihood(end);
            if diag.Likelihood(end) > bestlike
                bestlike = diag.Likelihood(end);
                bestgam = gam;
            end
        end
    end
    gammap = bestgam;
else
    likelihoods = [];
end

% kernelize the data and add bias
ktrials = utl_kernelize(trials,basis,kernel,gammap,polydegree);
ktrials = quickif(bias,[ones(size(ktrials,1),1) ktrials],ktrials);
% run the RVM
[param, hyperparam, diag] = hlp_diskcache('predictivemodels',@SparseBayes,likelihood,ktrials,targets,SB2_UserOptions(args1{:}),SB2_ParameterSettings(args2{:}));

% preselect relevant basis vectors
if ~strcmp(kernel,'linear')
    if bias
        feature_sel = param.Relevant-1;
        if feature_sel(1) == 0
            % bias was relevant, remove it from the basis vector selection (it will be added after kernelization)
            feature_sel = feature_sel(2:end);
        else
            % bias was not relevant, forget about it
            bias = 0;
        end
    else
        feature_sel = param.Relevant;
    end
    basis = basis(feature_sel,:);
else
    feature_sel = param.Relevant;
end

model = struct('sc_info',{sc_info},'classes',{classes},'basis',{basis},'param',{param},'hyperparam',{hyperparam},...
               'diag',{diag},'ptype',{ptype},'feature_sel',{feature_sel},'bias',{bias}, ...
               'kernel',{kernel},'gamma',{gammap},'gamma_likelihoods',likelihoods,'degree',{polydegree});
