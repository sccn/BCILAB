function model = ml_trainvote(trials, targets, votingscheme, learner, predictor, varargin)
% Internal meta-algorithm for voting. Used by other machine learning functions.
% Model = ml_trainvote(Trials, Targets, VotingScheme, BinaryLearner, BinaryPredictor, LearnerArguments...)
%
% Voting is a simple meta-algorithm that allows to use binary learners to learn a predictive model
% over multiple classes; Most ml_train* functions make implicit use of this function when the input
% has more than two classes. In the implemented 1-vs-1 voting schema, a binary predictive model is
% computed for every pair of classes, and during prediction, the probabilities that the predictors
% assign for class i are summed, for every i, and finally re-normalized.
%
% In:
%   Trials       : training data, as in ml_train
%
%   Targets      : target variable; should be a Nx1 array
%
%   VotingScheme : the voting scheme to use, can be one of:
%                   '1v1' : one-vs-one voting; if a binary classifier makes assumptions about the
%                           uniformity of its two classes (such as Linear Discriminant Analysis)
%                           and/or has trouble supporting unequal class priors, this is the
%                           preferred voter.
%
%                   '1vR' : one-vs-rest voting; if a binary classifier handles complex feature-space
%                           distributions well (i.e., is very robust) and produces reasonable
%                           probabilities (that are informed by potentially unequal class priors),
%                           then this is the preferred voter
%
%   BinaryLearner : binary learning function, e.g. one of the ml_train* functions
%                   accepts Trials & Targets as its first two arguments and optionally further
%                   arguments which can be supplied as Learner-Arguments
%
%   BinaryPredictor : binary prediction function, e.g. the corresponding ml_predict* function
%                     accepts Trials as its first argument and a model as second argument.
%
%   LearnerArguments : list of further arguments to be passed through to the learner
%                      (the first 2 are implicit, and a subset of Trials and Targets)
%
% Out:
%   Model   : a voted model
%             .classes is the target values per class
%             .voted is a cell array of binary models
%             .opts is the options
%
% Examples:
%   % learn one-vs-one voting arrangement of LDA classifiers, supporting multiple classes
%   model = ml_trainvote(trials,labels,'1v1',@ml_trainlda,@ml_predictlda) 
%
%   % as before, but learn a 1-vs-rest voting arrangement
%   model = ml_trainvote(trials,labels,'1vR',@ml_trainlda,@ml_predictlda)
%
%   % as before, but this time pass some options to the learning function
%   model = ml_trainvote(trials,labels,'1vR',@ml_trainlda,@ml_predictlda,[],'weight_bias',1,'weight_cov',1)
%
% See also:
%   ml_predictvote
%
%                           Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                           2010-06-25

% first get the classes for training
if iscell(trials) && iscell(targets)
    model.classes = unique([targets{:}]);
else
    model.classes = unique(targets);
end
    
switch votingscheme
    case '1v1'
        if iscell(targets)
            error('One-vs-one voting not yet support for multi-task learning problems.'); end
        for i=1:length(model.classes)
            for j=i+1:length(model.classes)
                % learn an i-vs-j model
                % get trial subset
                subset = targets == model.classes(i) | targets == model.classes(j);
                % learn restricted model
                if isnumeric(trials)
                    switch ndims(trials)
                        case 2
                            model.voted{i,j} = learner(varargin{:},'trials',trials(subset,:),'targets',targets(subset));
                        case 3
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,subset),'targets',targets(subset));
                        case 4
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,subset),'targets',targets(subset));
                        case 5
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,subset),'targets',targets(subset));
                        case 6
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,:,subset),'targets',targets(subset));
                        case 7
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,:,:,subset),'targets',targets(subset));
                        case 8
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,:,:,:,subset),'targets',targets(subset));
                        case 9
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,:,:,:,:,subset),'targets',targets(subset));
                        case 10
                            model.voted{i,j} = learner(varargin{:},'trials',trials(:,:,:,:,:,:,:,:,:,subset),'targets',targets(subset));
                        otherwise
                            error('At most 10-way feature tensors are supported by ml_trainvote.');
                    end
                else
                    % custom data
                    model.voted{i,j} = learner(varargin{:},'trials',utl_default_partitioner(trials,subset),'targets',targets(subset));
                end
            end
        end
    case '1vR'
        for i=1:length(model.classes)
            % learn model with i targets relabeled to 1 and remaining ones relabeled to 2
            if iscell(trials) && iscell(targets)
                model.voted{i} = learner(varargin{:},'trials',trials,'targets',cellfun(@(t)1+(t ~= model.classes(i)),targets,'UniformOutput',false));
            else
                model.voted{i} = learner(varargin{:},'trials',trials,'targets',1 + (targets ~= model.classes(i)));
            end
        end
    otherwise
        error(['unsupported voting scheme: ' votingscheme]);
end
   
if iscell(trials) && iscell(targets)
    if isnumeric(trials{1})
        model.feature_dims = ndims(trials{1})-1; end
else
    if isnumeric(trials)
        model.feature_dims = ndims(trials)-1; end
end

% and remember the options
model.opts = [{learner, predictor} varargin];
model.scheme = votingscheme;
