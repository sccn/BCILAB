function model = ml_trainhdca(varargin)
% Learn a linear predictive model by (regularized) Hierarchical Discriminant Component Analysis.
% Model = ml_trainhdca(Trials, Targets, Lambda, Options...)
%
% HDCA is a 2/3-level hierarchical classifier. The basic idea is that features are partitioned into
% "blocks" of a given size, and a classifier is trained for each block, followed by a classifier
% that acts on the linearly mapped output of the per-block classifiers. Optionally, this can be done
% for a number of given ranges of "channels" in Trials, and then fused with a third-level classifier.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable, as in ml_train
%
%   Lambda       : within-block regularization parameter, reasonable range: 0:0.1:1, greater is stronger
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
%   % learn a standard shrinkage HDCA model
%   model = ml_trainhdca(trials,targets);
%
%   % take unequal class priors into account for both the bias and the covariance matrix
%   model = ml_trainhdca(trials,targets,[],'weight_bias',1,'weight_cov',1);
%
%   % use a different type of regularization, which controls feature independence and requires cross-validation
%   model = utl_searchmodel({trials,target},'args',{{'lda',search(0:0.1:1),'regularization','independence'}})
%
%
% See also:
%   ml_predicthdca
%
% References:
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-04-03
        
args = arg_define([0 2],varargin, ...
    arg_norep('trials'), ...
    arg_norep('targets'), ...
    arg({'within_lambda','WithinLambda'}, [], [0 1], 'Optional regularization parameter. This is for the within-block classifiers. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence" (not necessary in "auto" mode).'), ...
    arg({'across_lambda','AcrossLambda'}, [], [0 1], 'Optional regularization parameter. This is for the across-block classifier. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence" (not necessary in "auto" mode).'), ...
    arg({'multimodal_lambda','MultimodalLambda'}, [], [0 1], 'Optional regularization parameter. This is for the multi-modal classifier. Reasonable range: 0:0.1:1 - greater is stronger. Requires that the regularization mode is set to either "shrinkage" or "independence" (not necessary in "auto" mode).'), ...
    arg({'regularization','Regularizer','Regularization'}, 'auto', {'none','auto','shrinkage','independence'}, 'Type of regularization. Regularizes the robustness / flexibility of covariance estimates. Auto is analytical covariance shrinkage, shrinkage is shrinkage as selected via plambda, and independence is feature independence, also selected via plambda.'), ...
    arg({'modality_ranges','ModalityRanges'},{},[],'Channel ranges for each modality. If empty, only one modality is assumed.','type','expression'), ...
    arg({'weight_bias','WeightedBias'}, false, [], 'Account for class priors in bias. If you do have unequal probabilities for the different classes, this should be enabled.'), ...
    arg({'weight_cov','WeightedCov'}, false, [], 'Account for class priors in covariance. If you do have unequal probabilities for the different classes, it makes sense to enable this.'), ...
    arg({'votingScheme','VotingScheme'},'1vR',{'1v1','1vR'},'Voting scheme. If multi-class classification is used, this determine how binary classifiers are arranged to solve the multi-class problem. 1v1 gets slow for large numbers of classes (as all pairs are tested), but can be more accurate than 1vR.'));

arg_toworkspace(args);

% find the class labels
classes = unique(targets);
if length(classes) > 2
    % learn a voting arrangement of models...
    model = ml_trainvote(trials, targets, votingScheme, @ml_trainhdca, @ml_predicthdca, varargin{:},'weight_bias',true);    %#ok<*NODEF>
elseif length(classes) == 1
    error('BCILAB:only_one_class','Your training data set has no trials for one of your classes; you need at least two classes to train a classifier.\n\nThe most likely reasons are that one of your target markers does not occur in the data, or that all your trials of a particular class are concentrated in a single short segment of your data (10 or 20 percent). The latter would be a problem with the experiment design.');
else
    % determine block sizes if necessary
    if ndims(trials) == 3
        if isempty(modality_ranges)
            modality_ranges = {1:size(trials,1)}; end        
    else
        error('This classifier requires 3d tensor features.');
    end
            
    % we call vanilla LDA to do the regressions - determine parameters
    lda_parameters = hlp_struct2varargin(args,'suppress',{'trials','targets','block_sizes','within_lambda','across_lambda','multimodal_lambda','modality_ranges'});
    
    % for each modality range
    for m=1:length(modality_ranges)
        range = modality_ranges{m};
        % for each block...
        for b=size(trials,2):-1:1
            blocktrials = reshape(trials(range,b,:),[],size(trials,3))';
            blockmodels{m}{b} = ml_trainlda(lda_parameters{:},'trials',blocktrials, 'targets',targets, 'lambda',within_lambda);
            blockpredictions{b} = blocktrials*blockmodels{m}{b}.w' - blockmodels{m}{b}.b';
        end
        layertrials = cat(2,blockpredictions{:});
        rangemodels{m} = ml_trainlda(lda_parameters{:}, 'trials',layertrials, 'targets',targets, 'lambda',across_lambda); %#ok<AGROW>
        rangepredictions{m} = layertrials*rangemodels{m}.w' - rangemodels{m}.b'; %#ok<AGROW>
    end
    toptrials = cat(2,rangepredictions{:});
    topmodel = ml_trainlda(lda_parameters{:}, 'trials',toptrials, 'targets',targets, 'lambda',multimodal_lambda);
    
    model = struct('modality_ranges',{modality_ranges}, 'blockmodels',{blockmodels}, 'rangemodels',{rangemodels}, 'topmodel',{topmodel});
end
