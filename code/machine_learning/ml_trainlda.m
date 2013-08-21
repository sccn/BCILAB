function model = ml_trainlda(varargin)
% Learn a linear predictive model by (regularized) Linear Discriminant Analysis.
% Model = ml_trainlda(Trials, Targets, Lambda, Options...)
%
% LDA is one of the simplest and oldest learning algorithms, originally introduced by Fisher in [1].
% Its basic assumption is that the data originates from two classes (extended to more via 1-vs-1
% voting), where the data in each class is distributed in the feature space according to a normal
% distribution. The shape of the distribution is assumed to be identical for both classes, and the
% relative weight/prior probability of each class is (by default) assumed to be identical, as well.
% An example would be data with features generated as a weighted sum of some latent variables which
% take on different values between conditions (but identical values within each condition), and
% which are superimposed with gaussian noise (which is identically distributed across conditions,
% e.g. from a large sum of independent random variables).
%
% Despite its simplicity, LDA assumes more structure in the data than is usually necessary (namely a
% certain normal distribution per class), and sometimes, more than what can be satisfactorily
% learned from the data, so that, even when the assumptions are fulfilled, the method is not
% guaranteed to give the optimal result. The estimation of the per-class covariance matrix is a
% notoriously data-hungry (and especially outlier-sensitive) step, and the main weakness of standard
% LDA.
%
% ml_trainlda offers three advanced covariance estimators, with different trade-offs to mitigate
% these problems. The variants 'shrinkage' and 'independence' each introduce a regularization
% parameter [2] which controls the complexity of the estimated matrix, and which need to be selected
% in a parameter search (which is orders of magnitude more time-consuming). The 'auto' variant
% computes the degree of regularization analyically in closed form, and is therefore both fast and
% (in some sense) optimal (but more restricted than 'independence') [3]. When enough trials are
% available, full covariance matrices are learned, but the less trials are given, the more the
% covariance estimates degrade to spherical (though well-formed) ones. Auto is the default for lda.
%
% While regularization can automatically control the complexity of a classifier, it is not a panacea
% which allows to add arbitrarily many features, since with each additional feature, the amount of
% structure (here:correlation) that can be captured in the remaining ones gets reduced.
%
% Within the toolbox, LDA is one of the bread-and-butter classifiers, and is worth trying in every
% reasonably simple classification task, especially for its speed. The major problem with LDA
% compared to other available classifiers is that it easier to break it with outliers than, for
% example, support vector machines or logistic regression. Another weakness is that the outputs are
% relatively primitive probability estimates, in contrast to, for example, logistic regression
% (somewhat better) or relevance vector machines (clearly better).
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Lambda       : optional regularization parameter, reasonable range: 0:0.1:1, greater is stronger
%                  requires that the regularization mode is set to either 'shrinkage' or 'independence' (default: [])
%           
%   Options  : optional name-value parameters to control the training details:
%              'regularization' -> 'shrinkage': covariance shrinkage, depends on plambda 
%                                  'independence': feature independence, depends on plambda
%                                  'auto': analytical covariance shrinkage, plambda is ignored (default)
%              'weight_bias' -> 0/1, take unequal class priors into account for bias calculation
%                               default: 0
%              'weight_cov' -> 0/1, take unequal class priors into account for covariance calculation
%                              default: 0
% Out:
%   Model   : a linear model; w is the linear weights, b is the bias; classes indicates the class labels which the model predicts
%
% Examples:
%   % learn a standard shrinkage LDA model
%   model = ml_trainlda(trials,targets);
%
%   % take unequal class priors into account for both the bias and the covariance matrix
%   model = ml_trainlda(trials,targets,[],'weight_bias',1,'weight_cov',1);
%
%   % use a different type of regularization, which controls feature independence and requires cross-validation
%   model = utl_searchmodel({trials,target},'args',{{'lda',search(0:0.1:1),'regularization','independence'}})
%
%
% See also:
%   ml_predictlda
%
% References:
%   [1] Fisher, R. "The use of multiple measurements in taxonomic problems."
%       Annals Eugen. 7 (1936), 188, 179.
%   [2] Friedman, J. "Regularized discriminant analysis." 
%       Journal of the American Statistical Association 84, 405 (1989), 175, 165.
%   [3] O. Ledoit and M. Wolf, "A well-conditioned estimator for large-dimensional covariance matrices"
%       J Multivar Anal, 88(2): 365-411, 2004.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03
        
arg_define([0 3],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'plambda','Lambda','lambda'}, [], [], 'Optional regularization parameter. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence" (not necessary in "auto" mode).'), ...
    arg({'regularization','Regularizer'}, 'auto', {'none','auto','shrinkage','independence'}, 'Type of regularization. Regularizes the robustness / flexibility of covariance estimates. Auto is analytical covariance shrinkage, shrinkage is shrinkage as selected via plambda, and independence is feature independence, also selected via plambda.'), ...
    arg({'weight_bias','WeightedBias'}, false, [], 'Account for class priors in bias. If you do have unequal probabilities for the different classes, this should be enabled.'), ...
    arg({'weight_cov','WeightedCov'}, false, [], 'Account for class priors in covariance. If you do have unequal probabilities for the different classes, it makes sense to enable this.'));

% find the class labels
classes = unique(targets);
if length(classes) > 2
    % learn a voting arrangement of models...
    model = ml_trainvote(trials, targets, '1vR', @ml_trainlda, @ml_predictlda, varargin{:},'weight_bias',true);    %#ok<*NODEF>
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    % pre-prune degenerate features
    retain = true(1,size(trials,2));
    for c = 1:length(classes)
        X = trials(targets==classes(c),:);
        n{c} = size(X,1);
        mu{c} = mean(X,1);
        v{c} = var(X,[],1);
        retain = retain & isfinite(mu{c}) & isfinite(v{c}) & (n{c}==1 | v{c} > eps);
    end    
    % apply feature mask...
    trials = trials(:,retain);
    % estimate distribution of each class...
    for c = 1:length(classes)
        X = trials(targets==classes(c),:);
        n{c} = size(X,1);
        mu{c} = mean(X,1);
        if n{c} == 1
            sig{c} = zeros(size(X,2));
        elseif strcmp(regularization,'auto')
            sig{c} = cov_shrink(X);
        else
            sig{c} = cov(X);
            if ~isempty(plambda) && ~strcmp(regularization,'none')
                % plambda-dependent regularization
                if strcmp(regularization,'independence')
                    sig{c} = (1-plambda)*sig{c} + plambda * diag(diag(sig{c}));
                elseif strcmp(regularization,'shrinkage')
                    sig{c} = (1-plambda)*sig{c} + plambda*eye(length(sig{c}))*abs(mean(diag(sig{c})));
                else
                    error('unknown regularization mode');
                end
            end
        end
    end
    
    ns = fastif(weight_cov,n,{1 1});
    nb = fastif(weight_bias,n,{1 1});
    % do the math
    mu_both = (mu{1}*nb{2} + mu{2}*nb{1}) / (nb{1}+nb{2});    
    sig_both = (sig{1}*ns{1} + sig{2}*ns{2}) / (ns{1}+ns{2});
    w = (mu{2} - mu{1}) / sig_both;
    w = w / (mu{2}*w' - mu_both*w');
    model = struct('w',{w}, 'b',{mu_both*w'}, 'classes',{classes},'featuremask',{retain});
end
