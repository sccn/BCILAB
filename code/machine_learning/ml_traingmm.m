function model = ml_traingmm(varargin)
% Learn a probabilistic non-linear predictive model using Gaussian Mixture Models.
% Model = ml_traingmm(Trials, Targets, Clusters, Options...)
% 
% Gaussian Mixture Models [1] are well-known generative statistical models that can be used for
% clustering, classification, and density estimation, whereas the provided implementation is
% primarily targeted at classification. The basics assumption of Gaussian mixture models is that the
% data (in each class) is a mixture of (relatively few) approximately gaussian-distributed distinct
% clusters. The optimization methods used for finding these clusters are usually very fragile
% (compared to, e.g., SVMs, where this problem is bypassed) and give reasonable results when the
% data is in fact a collection of distinct clusters, each of which deviates not too much from being
% gaussian, and when there are sufficient samples per cluster, w.r.t. to the total dimensionality of
% the data. These assumptions are usually only fulfilled when the data of one (or more) class(es)
% was derived from multiple distinct situations (e.g. walking, driving, waiting) and when either a
% large number of trials is available, or the dimensionality of the feature space is relatively low,
% i.e. the features come from a an effective adaptive feature-extraction stage.
%
% Several variants are supplied, among others "conventional" approaches using, e.g., expectation
% maximization, which require that the number of clusters is approximately known in advance. These
% are implemented via the GMMBAYES library [2]. Further, there are "advanced" implementations, where
% the number of clusters is not needed in advance, but learned from the data, usually using a
% variant of the Dirichlet process prior, which are implemented using the VDPGM [3] library.
%
% A benefit of GMMs is that the produced probability are of the highest quality among all methods
% implemented by the toolbox. For example, data points which are not in any of the previously
% learned classes typically receive equal probabilities under all classes (which cannot be expected
% from logistic regression or relevance vector machines). However, due to their restrictions,
% gaussian mixture models should be viewed as rather specialized classifiers, which should only be
% relied on when the assumptions are fulfilled.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Clusters     : expected/maximum number of clusters (required for some variants)
%                  default: 3
%
%   Options  : optional name-value parameters to control the training details:
%              'variant': one of several variants, including
%                        "nonparametric" methods:
%                        'avdp': accelerated variational Dirichlet process, ignores Clusters
%                        'vdp' : variational Dirichlet process iterative, using a per-weight prior, 
%                                   more or less uninformative, but numerically more robust than the 'vb' variant
%                        'bj'  : Blei & Jordan
%                        'cdp' : collapsed Dirichlet prior
%                        'csp' : collapsed stick-breaking
%                        conventional parametric methods:
%                        'vb'  : variational Bayes (exact # clusters)
%                        'em'  : expectation-maximization (exact # cluster)
%                        'fj'  : Figueiredo-Jain (max. # clusters)
%                        'gem' : greedy expectation-maximization (max. # clusters)
%              
%              'scaling': pre-scaling of the data (see hlp_findscaling for options) (default: 'std') 
%
% Out:
%   Model   : a multi-class gaussian mixture model; 
%             classes indicates the class labels which the model predicts
%
% Examples:
%   % learn a simple Gaussian Mixture model classifier (using an advanced method which learns the 
%   % number of clusters from the data)
%   model = ml_traingmm(trials,targets)
%
%   % use the variational Dirichlet process Gaussian mixture model classifier
%   model = ml_traingmm(trials,targets,[],'variant','vdp')
%
%   % use the simple Expectation-Maximization method, and assume 3 clusters of data points per class
%   model = ml_traingmm(trials,targets,3,'variant','em')
%
%   % like before, but search over possible values for the number of clusters
%   model = utl_searchmodel({trials,targets},'args',{{'gmm',search(1:5),'variant','em'}})
%
%   
% See also:
%   ml_predictgmm, gmmb_create
%
% References:
%   [1] Trevor Hastie and Robert Tibshirani, "Discriminant Analysis by Gaussian Mixtures"
%       Journal of the Royal Statistical Society. Series B (Methodological), Vol. 58, No. 1 (1996), pp. 155..176
%   [2] Nico Vlassis and A. Likas, "A greedy EM algorithm for gaussian mixture learning."
%       Neural Process. Lett. 15, 1 (2002), 77?87.
%   [3] Kenichi Kurihara, Max Welling and Nikos Vlassis "Accelerated Variational Dirichlet Mixture Models"
%       Advances in Neural Information Processing Systems 19 (NIPS 2006).
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-05

arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'clusters','NumClusters'}, 3, uint32([0 1000]), 'Expected/maximum number of clusters. Required for some variants (vb,em,fj,gem).'), ...
    arg({'variant','Variant'}, 'Accelerated Variational Dirichlet Process', {'Accelerated Variational Dirichlet Process','Variational Dirichlet Process','Blei-Jordan','Collapsed Dirichlet Process','Collapsed Stick-Breaking','Variational Bayes','Expectation-Maximization','Figuriedo-Jain','Greedy Expectation=Maximization'}, ['Variant to use. Non-parametric methods: Accelerated Variational Dirichlet Process, Variational Dirichlet Process, Blei & Jordan''s method, Collapsed Dirichlet Prior, ' ...
                                                                                              'Collasped Stick-Breaking construction. Conventional Methods: Variational Bayes, Expectation-Maximiation, Figuriedo-Jain, Greedy Expectation-Maximization.']), ...
    arg({'scaling','Scaling'}, 'std', {'none','center','std','minmax','whiten'}, 'Pre-scaling of the data. For the regulariation to work best, the features should either be naturally scaled well, or be artificially scaled.'), ...
    arg({'is_verbose','Verbose','verbose'},false,[],'Verbose outputs.'));

% scale the data
sc_info = hlp_findscaling(trials,scaling); %#ok<*NODEF>
trials = hlp_applyscaling(trials,sc_info);

% identify and remap the classes
classes = unique(targets);
% remap target labels to 1..k (bias is added later)
targ = targets;
for c=1:length(classes)
    targets(targ==classes(c)) = c; end
if length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.'); end

variant = hlp_rewrite(variant,'Accelerated Variational Dirichlet Process','avdp','Variational Dirichlet Process','vdp','Blei-Jordan','bj','Collapsed Dirichlet Process','cdp','Collapsed Stick-Breaking','csb','Variational Bayes','vb','Expectation-Maximization','em','Figuriedo-Jain','fj','Greedy Expectation=Maximization','gem');


% back up the rand state
if any(strcmp(variant,{'em','gem','fj'}))
    oldstate = rand('state'); %#ok<RAND>
    c = onCleanup(@() rand('state',oldstate)); %#ok<RAND>
    rand('seed',1337); %#ok<RAND>
end

switch variant
    case {'avdp','vdp','bj','bjrnd','cdp','csb','vb'}
        warning off MATLAB:divideByZero
        if any(strcmp(variant,{'avdp','vdp'}))
            opts = hlp_diskcache('predictivemodels',@feval,['mkopts_' variant]);
        else
            opts = hlp_diskcache('predictivemodels',@feval,['mkopts_' variant], clusters);
        end
        opts.get_q_of_z = 1;
        opts.seed = 1337;
        for c=1:length(classes)
            model.data{c} = trials(targets==c,:);
            if is_verbose
                model.class{c} = hlp_diskcache('predictivemodels',@vdpgm,model.data{c}',opts);
            else
                [t,model.class{c}] = evalc('hlp_diskcache(''predictivemodels'',@vdpgm,model.data{c}'',opts)'); %#ok<ASGLU>
            end
        end
    case 'em'
        warning off gmmb_em:data_amount
        model.class = hlp_diskcache('predictivemodels',@gmmb_create,trials, targets, 'EM', 'components', clusters);
    case 'fj'
        warning off gmmb_fj:data_amount
        model.class = hlp_diskcache('predictivemodels',@gmmb_create,trials, targets, 'FJ', 'Cmax', clusters);
    case 'gem'
        warning off gmmb_gem:data_amount
        model.class = hlp_diskcache('predictivemodels',@gmmb_create,trials, targets, 'GEM', 'Cmax', clusters);
    otherwise
        error('unknown variant specified');
end

% recover the rand state
if exist('oldstate','var')
    rand('state',oldstate); end %#ok<RAND>

model.classes = classes;
model.sc_info = sc_info;
model.variant = variant;
