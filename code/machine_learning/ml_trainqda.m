function model = ml_trainqda(varargin)
% Learn a non-linear predictive model by (regularized) Quadratic Discriminant Analysis.
% Model = ml_trainqda(Trials, Targets, Lambda, Kappa Options...)
%
% Quadratic discriminant analysis [1] is a (quadratic) generalization of linear discriminant
% analysis, and practically all remarks from ml_trainlda apply to this classifier, as well. However,
% QDA is siginificantly more prone to overfitting than is LDA, since it relaxes the assumption that
% the covariance structure of each class is identical. For this reason, it has twice as many
% parameters, and is also approx. twice as slow. Aside from that, QDA is the simplest non-linear
% classifier and is interesting especially when it is known that one distribution is shaped
% significantly different from the other (e.g. in one-vs-rest votings, or in baseline-versus-signal
% classification).
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Lambda       : optional per-class covariance regularization parameter, reasonable range: 0:0.1:1, greater is stronger
%                  requires that the regularization option is set to either 'shrinkage' or 'independence' (default: [])
%   Kappa        : between-class covariance regularization parameter, reasonable range: 0:0.1:1, greater is stronger (default: [])
%
%   Options  : optional name-value parameters to control the training details:
%              'weight_bias' -> 0/1, take unequal class priors into account for bias calculation
%                               default: 0 -- note: this parameter is currently ignored
%              'weight_cov' -> 0/1, take unequal class priors into account for covariance calculation
%                              default: 0
%              'regularization' -> 'shrinkage': covariance shrinkage, depends on lambda 
%                                  'independence': feature independence, depends on lambda
%                                  'auto': analytical covariance shrinkage, lambda is ignored (default)
% Out:
%   Model   : a quadratic model; q is the quadratic weights, l is the linear weights, b is the bias; 
%             classes indicates the class labels which the model predicts; sc_info holds scale parameters
%
% Examples:
%   % learn a shrinkage QDA model
%   model = ml_trainqda(trials,targets);
%
%   % take unequal class priors into account for both the bias and the covariance matrix
%   model = ml_trainqda(trials,targets,[],[],'weight_bias',1,'weight_cov',1);
%
%   % learn a shrinkage QDA model which blends covariance estimates between classes, using a 
%   % regularization parameter
%   model = utl_searchmodel({trials,target},'args',{{'qda',[],search(0:0.1:1)}})
%
%   % use a different type of regularization, which controls feature independence
%   model = utl_searchmodel({trials,target},'args',{{'qda',search(0:0.1:1),[],'regularization','independence'}})
%   
%   % as before, but search for the optimal blending parameter at the same time (slow!)
%   model = utl_searchmodel({trials,target},'args',{{'qda',search(0:0.1:1),search(0:0.1:1),'regularization','independence'}})
%
% See also:
%   ml_predictqda
%
% References:
%   [1] Friedman, J. "Regularized discriminant analysis." 
%       Journal of the American Statistical Association 84, 405 (1989), 175, 165.
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03
    
arg_define([0 4],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'lambda','Lambda'}, [], [], 'Within-class covariance regularization parameter. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence".'), ...
    arg({'kappav','Kappa','kappa'}, [], [], 'Between-class covariance regularization parameter. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence".'), ...
    arg({'regularization','Regularizer'}, 'auto', {'auto','shrinkage','independence'}, 'Regularization type. Regularizes the robustness / flexibility of covariance estimates. Auto is analytical covariance shrinkage, shrinkage is shrinkage as selected via lambda, and independence is feature independence, also selected via lambda.'), ...
    arg({'weight_bias','WeightedBias'}, false, [], 'Account for class priors in bias. If you do have unequal probabilities for the different classes, this should be enabled.'), ...
    arg({'weight_cov','WeightedCov'}, false, [], 'Account for class priors in covariance. If you do have unequal probabilities for the different classes, it makes sense to enable this.'));


% find the class labels
classes = unique(targets);
if length(classes) > 2
    % learn a voting arrangement
    model = ml_trainvote(trials,targets,'1vR',@ml_trainqda,@ml_predictqda,varargin{:}); %#ok<*NODEF>
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    % pre-prune degenerate features
    retain = true(1,size(trials,2));
    for c = 1:length(classes)
        X = trials(targets==classes(c),:);
        n{c} = size(X,1);
        mu{c} = mean(X);
        v{c} = var(X);
        retain = retain & isfinite(mu{c}) & isfinite(v{c}) & v{c} > eps;
    end    
    % apply feature mask...
    trials = trials(:,retain);
    
    % scale the data (to keep the numerics happy)
    sc_info = hlp_findscaling(trials,'std');
    trials = hlp_applyscaling(trials,sc_info);
        
    % estimate distribution of each class...
    for c = 1:length(classes)
        % get mean and covariance
        X = trials(targets==classes(c),:);
        n{c} = size(X,1);
        mu{c} = mean(X);
        if strcmp(regularization,'auto')
            sig{c} = cov_shrink(X);
        else
            sig{c} = cov(X);
            if ~isempty(lambda)
                % lambda-dependent regularization
                if strcmp(regularization,'independence')
                    sig{c} = (1-lambda)*sig{c} + lambda * diag(diag(sig{c}));
                elseif strcmp(regularization,'shrinkage')
                    sig{c} = (1-lambda)*sig{c} + lambda*eye(length(sig{c}))*abs(mean(diag(sig{c})));
                else
                    error('unknown regularization mode');
                end
            end
        end
    end
    
    if ~isempty(kappav)
        % regularize the covariance matrices w.r.t. each other
        ns = fastif(weight_cov,n,{1 1});
        sigma = (sig{1}*ns{1} + sig{2}*ns{2}) / (ns{1}+ns{2});
        sig{1} = sig{1} * (1-kappav) + sigma*kappav;
        sig{2} = sig{2} * (1-kappav) + sigma*kappav;
    end
    
    % compute the model
    model = struct('c',{1/2*(logdet(sig{1})-logdet(sig{2})) + 1/2*((mu{1}/sig{1})*mu{1}' - (mu{2}/sig{2})*mu{2}')}, ...
        'l',{mu{1}/sig{1} - mu{2}/sig{2}}, 'q',{-1/2*(inv(sig{1}) - inv(sig{2}))}, 'sc_info',{sc_info}, 'classes',{classes}, 'featuremask',{retain});
end
