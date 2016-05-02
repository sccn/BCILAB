function model = ml_trainsvmperf(varargin)
% Learn a linear or non-linear predictive model by Support Vector Machines, using SVMperf.
% Model = ml_trainsvmperf(Trials, Targets, Cost, Options...)
%
% SVMperf [1] is a compehensive package providing linear and non-linear (kernel) Support Vector
% Machines [2] to do classification. Among the provided SVM implementations, it is the preferred
% method whenever more than just mis-classification rate is reqired as optimization criterion (more
% precisely: as loss function).
%
% Support Vector Machines are available for various machine learning tasks (SVMperf supports
% classification, SVMlight additionally supports regression and ranking), and are generally an
% excellent (fast and/or robust, depending on regularization choice) solution. Kernel SVMs are a
% special variant which allow for state-of-the-art non-linear learning (depending on the kernel
% type, the kernel parameter needs to be searched, in addition to the cost parameter). The major
% drawback of support vector machines is that they do not produce useful probability outputs, which
% restricts their use in decision-making setups. The other drawback is the need for regularization,
% which can be complicated and time-consuming (e.g. incurring nested cross-validation).
% 
% When the data is known to require a non-linear classifier (or regressor), Support Vector Machines,
% Relevance Vector Machines (ml_trainrvm) and Hierarchical Kernel Learning (ml_trainhkl) are
% currently the only methods provided in the toolbox, and their major difference is whether they
% give probabilistic outputs and how much regularization they require (and how well they are able to
% pick out relevant information in high-dimensional feature spaces).
%
% SVMperf supports a variety of kernels (rbf being the most useful one in most cases) a variety of
% advanced loss functions, and various computational methods/variants that primarily differ in speed
% (depending on data size, number of features, etc.) [3,4,5,6]. SVMperf is not free for
% non-commercial uses.
%
% In:
%   Trials   : training data, as in ml_train
%
%   Targets  : target variable, as in ml_train
%
%   Cost     : regularization parameter, reasonable range: 2.^(-5:2:15), greater is stronger
%
%   Options  : optional name-value parameters to control the training details
%               'loss': Loss function to use (0..5,10)
%                       * 0/'ZO'    Zero/one loss: 1 if vector of predictions contains error, 0 otherwise.
%                       * 1/'F1'    F1: 100 minus the F1-score in percent.
%                       * 2/'Err'   Errorrate: Percentage of errors in prediction vector.
%                       * 3/'PRBEP' Prec/Rec Breakeven: 100 minus PRBEP in percent.
%                       * 4/'Prec'  Prec@k: 100 minus precision at k in percent.
%                       * 5/'Rec'   Rec@k: 100 minus recall at k in percent.
%                       * 10/'ROC'  ROCArea: Percentage of swapped pos/neg pairs (i.e. 100 - ROCArea).
%               'method': choice of structural learning algorithm (0..9), see arg declaration below for details
%               'variant': CPSP variant for sparse kernel training (0,1,2,4, if method is 9), see arg declaration below for details
%               'slacknorm': L-norm to use for slack variables (1,2)
%               'loss_rescaling': Rescaling method to use for loss (1,2)
%               'sparse_basis': number of basis functions for sparse kernel approximation (e.g. 500)
%               'restarts': number of restarts during sparse kernel approximation (0+)
%              kernel parameters:
%               'kernel': ptype of kernel function (linear,poly,rbf,sigmoid,user); (default: 'rbf')
%               'gamma': parameter gamma in rbf kernel; reasonable search range: 2.^(-16:2:4) (default: 0.3)
%               'd': parameter d in polynomial kernel (default: 3)
%               's': parameter s in sigmoid/poly kernel (default: 1)
%               'r': parameter c in sigmoid/poly kernel (default: 1)
%               'u': parameter of user-defined kernel (default: '1')
%              misc options:
%               'eps': tolerance (e.g., 0.1)
%               'bias': bias present? (0,1, default:1)
%               'pr_k': k in Precision/Recall loss functions
%               'shrinking_heuristic': whether shrinking heuristic is used (0,1)
%               'quiet': quiet mode (0,1)
%               'scaling': pre-scaling, see hlp_findscaling (default: 'std')
%
% Out:
%   Model   : a linear model; 
%             classes indicates the class labels which the model predicts
%             sc_info is the scaling info
%
% Notes:
%   uses the SVM-struct learning module: SVM-perf, V3.00, 15.07.2009
%     includes SVM-struct V3.10 for learning complex outputs, 14.08.08
%     includes SVM-light V6.20 quadratic optimizer, 14.08.08
%     SVM-perf, SVM-struct, SVM-light copyright (C) Thorsten Joachims, thorsten@joachims.org
%
% Examples:
%   % learn a quick and dirty SVM model (without parameter search)
%   model = ml_trainsvmperf(trials,labels)
%
%   % learn an SVM model by searching over the cost parameter (note: quite slow)
%   model = utl_searchmodel({trials,labels},'args',{{'svmperf',search(2.^(-5:2:15))}})
%
%   % as before, but also search over the kernel scale parameter (note: really slow)
%   model = utl_searchmodel({trials,labels},'args',{{'svmperf',search(2.^(-5:2:15)),'gamma',search(2.^(-16:2:4))}})
%
%   % as before, but use F1 loss
%   model = utl_searchmodel({trials,labels},'args',{{'svmperf',search(2.^(-5:2:15)),'loss','F1','gamma',search(2.^(-16:2:4))}})
%
%   % as before, but use a linear kernel (no need to search over gamma, then)
%   model = utl_searchmodel({trials,labels},'args',{{'svmperf',search(2.^(-5:2:15)),'kernel','linear'}})
%
%
% See also:
%   ml_predictsvmperf
%
% References:
%   [1] T. Joachims, "A Support Vector Method for Multivariate Performance Measures", 
%       Proceedings of the International Conference on Machine Learning (ICML), 2005.
%   [2] Schoelkopf, B., and Smola, A. "Learning with Kernels: Support Vector Machines, Regularization, Optimization, and Beyond"
%       (Adaptive Computation and Machine Learning). The MIT Press, Dec. 2001.
%   [3] T. Joachims, "Training Linear SVMs in Linear Time", 
%       Proceedings of the ACM Conference on Knowledge Discovery and Data Mining (KDD), 2006.
%   [4] T. Joachims, Chun-Nam Yu, "Sparse Kernel SVMs via Cutting-Plane Training",
%       Proceedings of the European Conference on Machine Learning (ECML), 2009.
%   [5] I. Tsochantaridis, T. Joachims, T. Hofmann, and Y. Altun, "Large Margin Methods for Structured and Interdependent Output Variables", 
%       Journal of Machine Learning Research (JMLR), Vol. 6(Sep):1453-1484, 2005.
%   [6] T. Joachims, T. Finley, Chun-Nam Yu, "Cutting-Plane Training of Structural SVMs", 
%       Machine Learning Journal, to appear.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-04
dp;

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'cost','Cost'}, search(2.^(-5:2:15)), [0 2^-7 2^15 Inf], 'Regularization parameter. Reasonable range: 2.^(-5:2:15), greater is stronger. By default, it is the (average of x*x)^-1.','cat','Core Parameters'), ...
    arg({'loss_function','Loss','loss'}, 'Err', {'Z0','F1','Err','PRBEP','Prec','Rec','ROC'}, 'Loss function to use. Zero/one loss, F1-score, error rate, Precision/Recall at break-even point, Precision at k, Recall at k, ROC area (percentage of swapped pos/neg pairs). Note that some losses x are actually 1-x (since their natural formulation is "positive" rather than in terms of loss).','cat','Core Parameters'), ...
    arg({'method','LearningAlgo'}, 'Custom', {'n-Slack','n-SlackShrink','1-SlackPrimal','1-SlackDual','1-SlackDualCCache','Custom'}, 'Structural learning algorithm. N-slack algorithm, n-slack algorithm with shrinking heuristic, 1-slack algorithm (primal), 1-slack algorithm (dual), 1-slack algorithm (dual) with constraint cache, or the custom algorithm (see machine_learning/ml_trainsvmperf for references).','cat','Core Parameters'), ...
    arg({'variant','Method9CPSP'}, 'FixedPoint', {'off','SubsetHeurstic','FixedPoint','FixedPointSubsetHeuristic'}, 'Variant of the CPSP algorithm. For sparse kernel training; Either do not use CPSP, or CPSP using subset selection for preimages via 59/95 heuristic, or CPSP using fixed point search (RBF kernel only), or CPSP using fixed point search with starting point via 59/95 heuristic (RBF kernel only).','cat','Core Parameters'), ...
    arg({'kernel','Kernel'}, 'rbf', {'linear','rbf','poly','sigmoid','user'}, 'Kernel type. Linear, or Non-linear kernel types: Radial Basis Functions (general-purpose),  Polynomial (rarely preferred), Sigmoid (usually overly simple), User (user-defined kernel from kernel.h).','cat','Core Parameters'), ...
    arg({'g','RBFScale','gamma'}, search(2.^(-16:2:4)), [], 'Scaling parameter of the RBF kernel. Should match the size of structures in the data; A reasonable range is 2.^(-16:2:4).','cat','Core Parameters'), ...
    ...
    arg({'d','PolyDegree'}, 3, uint32([1 10]), 'Degree for the polynomial kernel.','cat','Miscellaneous'), ...
    arg({'s','SigmoidPolyScale'}, 1, [], 'Scale of sigmoid/polynomial kernel.','cat','Miscellaneous'), ...
    arg({'r','SigmoidPolyBias'}, 1, [], 'Bias of sigmoid / polynomial kernel.','cat','Miscellaneous'), ...
    arg({'u','UserParameter'}, '1', [], 'User-defined kernel parameter.','cat','Miscellaneous','type','char','shape','row'), ...
    arg({'pr_k','LossK'}, 0, [0 1], 'Fraction of positive examples k in Prec@k and Rec@k. Zero indicates to use 0.5x for Prec@k and 2x for Rec@k, with x being the number of positive examples in the training set.','cat','Miscellaneous'), ...
    arg({'slacknorm','SlackNorm'}, 'l1', {'l1','l2'}, 'L-norm to use for slack variables. Either L1-norm, or squared slacks.','cat','Miscellaneous','guru',true), ...
    arg({'loss_rescaling','LossRescaling'}, 'margin', {'slack','margin'}, 'Rescaling method to use for loss.','cat','Miscellaneous','guru',true), ...
    arg({'sparse_basis','NumBasisFunc'}, uint32(500), [0 Inf], 'Mumber of basis functions. For sparse kernel approximation.','cat','Miscellaneous','guru',true), ...
    arg({'restarts','NumRestarts'}, 0, uint32([0 1000]), 'Number of restarts. Number of times to recompute the sparse kernel approximation and restart the optimizer.','cat','Miscellaneous','guru',true), ...
    arg({'shrinking_heuristic','ShrinkingHeuristic'}, true, [], 'Use shrinking heuristic in custom aglorithm. Only for linear kernel and errorrate loss.','cat','Miscellaneous','guru',true), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.','cat','Miscellaneous'), ...
    arg({'epsi','Epsilon','eps'}, 0.1, [], 'Tolerated solution accuracy.','cat','Miscellaneous'), ...
    arg({'bias','Bias'}, false, [], 'Include a bias term. Only implemented for linear kernel.','cat','Miscellaneous'), ...
    arg({'votingScheme','VotingScheme'},'1vR',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'), ...
    arg({'verbose','Verbose'}, false, [], 'Show diagnostic output.','cat','Miscellaneous'));

if is_search(cost)
    cost = 1; end
if is_search(g)
    g = 0.3; end

% scale the data
sc_info = hlp_findscaling(trials,scaling); %#ok<*NODEF>
trials = hlp_applyscaling(trials,sc_info);

% find the class labels
classes = unique(targets);
if length(classes) > 2
    % in this case we use the voter
    model = ml_trainvote(trials,targets,votingScheme,@ml_trainsvmperf,@ml_predictsvmperf,varargin{:});
elseif length(classes) == 1
	error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    
    % remap target labels to -1,+1
    targets(targets==classes(1)) = -1;
    targets(targets==classes(2)) = +1;

    % rewrite sme string args to numbers
    loss_function = hlp_rewrite(loss_function,'ZO',0,'F1',1,'Err',2,'PRBEP',3,'Prec',4,'Rec',5,'ROC',10);
    method = hlp_rewrite(method,'n-Slack',0,'n-SlackShrink',1,'1-SlackPrimal',2,'1-SlackDual',3,'1-SlackDualCCache',4,'Custom',9);
    kernel = hlp_rewrite(kernel,'linear',0,'poly',1,'rbf',2,'sigmoid',3,'user',4);
    variant = hlp_rewrite(variant,'off',0,'SubsetHeurstic',1,'FixedPoint',2,'FixedPointSubsetHeuristic',4);
    slacknorm = hlp_rewrite(slacknorm,'l1',1,'l2',2);
    loss_rescaling = hlp_rewrite(loss_rescaling,'slack',1,'margin',2);

    % build the args
    args = sprintf('-c %f -v %d -y %d -p %d -o %d -l %d -w %d -e %f -t %d --p %d --k %d --r %d --s %d -d %d -g %f -s %f -r %f -u %s --b %d --i %d', ...
        cost,verbose,verbose,slacknorm,loss_rescaling,loss_function,method,epsi,kernel,pr_k,sparse_basis,restarts,shrinking_heuristic,d,g,s,r,u,bias,variant); % note: verbose deliberately shows up twice

    % run the command
    try
        model = hlp_diskcache('predictivemodels',@svmperflearn,trials,targets,args);
    catch e
        if strcmp(e.identifier,'MATLAB:UndefinedFunction')
            error('The SVMperf library has not been compiled for your platform. Please consider using one of the other SVM implementations instead.'); 
        else
            rethrow(e);
        end
    end
    model.sc_info = sc_info;
    model.classes = classes;
end
