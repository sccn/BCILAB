function model = ml_train(varargin)
% Learn a predictive model from the given data and parameters.
% Model = ml_train(Data, Arguments)
%
% This function dispatches, depending on the first cell in Arguments, to one of the
% machine_learning/ml_train* functions (e.g., ml_trainlda, etc.), and may perform some additional
% data formatting depending on inputs and outputs. It is to be preferred over calling the ml_train*
% functions directly (especially since some callees may rely on the additional work done by
% ml_train).
%
% The function is the general tool to be applied when a predictive model shall be learned from a
% (usually unstructured, but labeled) set of trials [1,3]. Each trial is typically fed to the
% function as a vector of previously extracted features (called a feature vector), as a generic
% format understood by a large variety of different methods [2]. Each trial typically has a label
% (or target value), which may be categorical (a number out of a small set of numbers), real-valued
% (or continuous: any number), multivariate (a vector of numbers), or structured (an arbitrary
% cell). If empty labels are supplied, the goal is usually not to predict any labels, but to learn
% the probablity distribution of the supplied data.
%
% The resulting predictive model (usually a MATLAB struct()) is termed, depending on the label type,
% a 'classifier', 'regressor', or 'density estimator'. It captures in its parameters the necessary
% information to be able to predict the desired target value from a new instance of data (assumed to
% originate from the same distribution as the training data). It is sometimes called the 'predictive
% function', and the term 'model' is also frequently used in the literature to refer to classes of
% such functions (with free parameters).
%
% Chosing the right method for a given task is a difficult problem, and requires correct assumptions
% about the structure of the data, and often further experimentation. The toolbox provides a
% relatively broad set of methods from different branches of machine learning, among others
% statistics (lda,qda,logreg), optimization (svm*,logreg,dal), support vector machines
% (svmlight,svmlinear,svmperf,rvm), (Bayesian) probability theory (gmm,logreg,rvm) and meta-learning
% (vote). Prominent fields currently not covered are neural networks, fuzzy inference, and
% evolutionary methods.
%
% The ml_train* functions are named according to the fundamental approach, but typically several
% variants are provided by each, which themselves may have multiple parameters and which may be
% implemented with different libraries. The overall approach can often be chosen based on
% assumptions about the data, whereas some parameters may have to be selected by parameter search
% (preferred for reasons explained in offline_analysis/bci_train), or are otherwise adapted. In the
% context of the toolbox, it is usually not enough to select the right learning function, but
% instead, the entire data pipeline, usually starting with signal (pre-)processing, followed by
% feature extraction and concluded by machine learning needs to be determined. Extensive fully
% customizable default setups ("paradigms") are provided at higher levels (paradigms/para_*,
% offline_analysis/bci_train).
%
% In:
%   Data    :   * either a cell array of {Trials} or {Trials, Targets} or {Trials, Targets,
%                 Weights}, with 
%                 - Trials a [NxF] numeric array (N... number of training instances,
%                   F... number of feature dimensions) or a [UxVxWx...xN] Nd array in special cases 
%                   (N... number of training instances, U,V,W,... feature tensor dimensions).
%                   If not otherwise possible, it is allowed to pass trials in any custom
%                   format, as long as it is being recognized by the machine learning function of
%                   interest.
%                 - Targets one of the following:
%                    * a [NxT] array (N... number of training instances, T... number of target 
%                      dimensions, T usually 1)
%                    * a {N} cell array (N... number of training instances, contents unspecified)
%                      ... in this case, the method selected in Options must support structured 
%                      prediction
%                    * an empty array, indicating that the probability density of the input data is 
%                      to be estimated
%                    * a cell array describing a set of probability distributions, either discrete
%                      or continuous (see ml_predict for examples)
%                 - Weights a [Nx1] numeric array of non-negative per-trial trial weights (N...
%                   number of training instances)
%               * or, alternatively, a cell array of {{[N1xF],[N2xF], [N3xF], ...}}, with
%                   Ni the number of instances for class i (possibly empty)
%                   F the number of feature dimensions
%
%   Arguments :   cell array of further arguments for learning the model; The first argument in
%                 Options is considered the name of the learning method, e.g. 'lda', 'logreg', and
%                 specifies the ml_train<name> learning function to be called all other options are
%                 passed unmodified into the learning function
% Out:
%   Model   : predictive model, can later be used in conjunction with ml_predict
%              Model.args is the arguments to ml_train
%              Model.model is the predictive model returned by the appropriate ml_train* function
%              where applicable,
%                 Model.model.b is the bias
%                 Model.model.w is the weights
%
% Examples:
%   % assuming a feature matrix for 160 trials a 10 dimensions, and a corresponding label vector
%   % size(trials) -> [160,10]; size(targets) -> [160,1]
%
%   % learn LDA predictor from some training data
%   model = ml_train({trials,targets},{'lda'})
%
%   % learn GMM predictor with some additional parameters
%   model = ml_train({trials,targets},{'gmm' 0.4 'xxx'})
%
%   % learn SVM predictor for three-class data in trials1/2/3
%   model = ml_train({{trials1,trials2,trials3}},{'svm'})
%
%   % learn GMM density estimator
%   model = ml_train({trials,[]}, {'gmm'})
%
%   % special: re-parse the list of supported plugin functions (affects GUIs displayed for ml_train)
%   ml_train('update')
%
%   % apply the learned model to some data
%   % predictions = ml_predict(model,trials)
%
% See also:
%   ml_predict
%
% Notes:
%   Weighted learning is only supported by a fraction of the ml_train* functions.
%
% References:
%   [1] Bishop, C. M. "Pattern Recognition and Machine Learning."
%       Information Science and Statistics. Springer, 2006.
%   [2] Hastie, T., Tibshirani, R., and Friedman, J. H. "The elements of statistical learning (2nd Ed.)."
%	    Springer, 2009.
%   [3] MacKay, D. J. C. "Information theory, inference, and learning algorithms."
%       Cambridge University Press, 2003.
%
% See also:
%   ml_traindal, ml_trainglm, ml_trainhkl, ml_trainlogreg, ml_trainrvm (built-in linear/logistic regression family),
%   ml_traingauss, ml_traingmm (built-in Gaussian mixture family),
%   ml_trainlda, ml_trainqda (built-in discriminant analysis family),
%   ml_trainsvm, ml_trainsvmlight, ml_trainsvmperf (built-in Support Vector
%   Machines family), ml_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-03
dp;

if ~isequal(varargin,{'update'})
    
    arg_define([0 2],varargin, ...
        arg_norep('data'), ...
        arg_subswitch({'learner','Learner'}, 'lda', list_learners(), 'Machine learning function. Applied to the data (features) produced within the paradigm; this is usually the last (and most adaptive) step in the processing from raw data to model or prediction.'));
    
    if ~iscell(data) %#ok<*USENS>
        error('Data must be a cell array.'); end
    
    if length(data) == 3
        % we have the {Trials,Targets,Weights} data format
        % this is passed on as ml_train*(trials,{targets,weights}, ...)
        trials = data{1};
        targets = data(2:3);
    elseif length(data) == 2
        % we have the {Trials,Targets} data format
        % this is passed on as ml_train*(trials,targets, ...)
        trials = data{1};
        targets = data{2};
    elseif length(data) == 1 
        if iscell(data{1}) && length(unique(cellfun(@(d)size(d,2),data{1}))) == 1
            % we have a pack of per-class cell arrays: concatenate & construct labels
            trials = cat(1,data{1}{:});
            for k=1:length(data{1})
                targets{k} = k*ones(size(data{1}{k},1),1); end
            targets = cat(1,targets{:});
        else
            trials = data{1};
            targets = [];
        end
    else
        % unknown data specification
        error('Data is specified in an unsupported format (must be a 1/2/3-element cell array).');
    end
    
    if isempty(targets)
        model.model = feval(['ml_train' learner.arg_selection],'trials',trials,learner);
    else
        if ~iscell(targets) && length(unique(targets)) == 1
            disp('WARNING: This training set contains only one class - the subsequent learning phase will likely fail.'); end
        model.model = feval(['ml_train' learner.arg_selection],'trials',trials,'targets',targets,learner);
    end
    model.args = learner;
    
else
    % update the list of learners
    list_learners(true);
end



function learners = list_learners(update_list)
dp;
% list all the learning functions in code/machine_learning/
global tracking;
persistent memo;
if isempty(memo) || exist('update_list','var') && update_list
    memo = {};
    ml_paths = {'functions:/machine_learning/ml_train*.m','home:/.bcilab/code/machine_learning/ml_train*.m'};
    if ~isempty(tracking.paths.private_path)
        ml_paths = [ml_paths {'private:/code/machine_learning/ml_train*.m'}]; end
    for p = ml_paths
        modules = dir(env_translatepath(p{1}));
        names = setdiff({modules.name},{'ml_train.m','ml_trainvote.m'});
        tags = cellfun(@(n) n(9:end-2),names,'UniformOutput',false);
        funcs = cellfun(@(n) str2func(n(1:end-2)),names,'UniformOutput',false);
        tmp = [tags; funcs];
        tmp = tmp(:)';
        memo = [memo tmp];
    end
end
learners = memo;
