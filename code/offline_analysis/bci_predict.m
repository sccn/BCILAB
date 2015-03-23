function [prediction, measure, stats, target] = bci_predict(varargin)
% Apply a predictive model to some data and optionally measure its performance.
% [Prediction, Loss, Statistics, Target] = bci_predict(Model,Data,TargetMarkers,EvaluationMetric)
%
% Let a model make predictions of cognitive state on a data set (usually with known target values),
% yielding one prediction per trial in the data, as well as the average empirical loss, if the data
% set has an attached target variable.
%
% bci_predict is used to obtain and/or evaluate the outputs of a predictive model on a recorded
% (usually raw) data set. Predictions are made for every trial in the data set, where trials are
% derived from the data in the same way as was used when the model was originally calibrated. For
% example, if the model was calibrated w.r.t. events of certain types in the calibration set (using
% the 'events' parameter of the paradigm), trials will be extracted for these same events in the
% supplied data set, as well. In the typical case where a target variable can be derived from the
% data set in this way (e.g. via the 'event' clause passed to bci_train - see bci_train), then the
% predictions of the model are compared with the defined target values, and an average loss is
% computed (and returned as the second output). This loss measure has the same meaning as in
% bci_train, and can be customized via the 'metric' parameter, either as one of the predefined
% losses (described in detail in machine_learning/ml_calcloss), or as a user-supplied custom
% function. By default, the most appropriate loss is automatically selected depending on the types
% of target and prediction values (mis-classification rate for categorical target values,
% mean-square error for continuous target values, etc.), so that user intervention is rarely needed.
%
% For advanced use, the actual predictions for every trial in the data are returned, as well as the
% defined target values in the data set (if available). The format of the targets depends on the
% procedure by which they were assigned to the data (e.g., when the 'events' clause was used to
% assign target values to certain events, targets is a numerical vector of class numbers, e.g. [1 2
% 2 1 4 1 3 2 3 2 1] for some hypothetical 4-class data set with 11 trials). The format of the
% predictions depends on the format of the targets (e.g., whether they were categorical values or
% continuous values, usually the former), and the capabilities of the 'learner' component of the
% paradigm. A detailed exposition of the possibilities is given in machine_learning/ml_predict.
%
% Almost all learners produce discrete probability distributions for categorical targets, so this is
% the format of the predictions in 90% of the situations. A sequence of discrete probability
% distributions (one per trial) as produced by bci_predict is formatted as a MATLAB cell array with
% the three entries {'disc', Probabilities, Classes}. 'disc' is the format tag of the predictions
% and indicates discrete probabilities, Probabilities is an array of probabilities for each trial
% and class, with # of trials rows and # of classes columns, i.e. it has size [NxC], and Classes is
% a column vector which specifies the desired target value for each of the C classes, in the same
% order as the rows in Probabilities. For this reason, the expected target value for any trial is
% pred{2}*pred{3}, and the most likely target value for any trial can be obtained as
% pred{3}(argmax(pred{2}')), if pred are the predictions return by bci_predict.
%
% In addition, further statistics are returned, which may contain additional values (e.g. false
% positive rate in the case of binary classification).
%
% Example:
%
%   % load a calibration data set
%   calib = io_loadset('data sets/john/gestures.eeg')
%   % train a model (see bci_train for an explanation of this stage)...
%   [loss,model,stats] = bci_train({'data',calib, 'paradigm',@para_csp}, 'events',{'left-imag','right-imag','rest'})
%   
%   % load a test data set
%   testdata = io_loadset('data sets/john/gestures2.eeg')
%
%   % and apply the model on that data set to get predictions and a loss estimate...
%   [predictions,loss2,stats2,targets] = bci_predict(model,testdata);
%
% The application of the model results in predictions for every of its trials, an empirical average
% loss measure, statistics, and the original target values for every trial in the test data. This is
% under the assumption that the file gestures2 contains markers 'left-imag', 'right-imag', 'rest',
% just like the first data set, so that a loss (deviation from desired outcomes) can be computed.
% Note that applying a model to its own calibration data set gives overly optimistic results, and is
% usually scientifically invalid.
%
%
% A less common use case is to exchange the default preprocessing that is applied by a paradigm by
% some alternative preprocessing (for example, one which slighty deviates from the model's settings,
% or a manually implemented one, or the by the preprocessing of an entirely different model). To
% this end, not the raw data is passed to bci_predict, but instead data that was manually processed,
% for example using bci_preproc (see bci_preproc). In the follwing example, a test data set is
% loaded (as in the previous sample), but it is manually preprocessed using the model's defined
% preprocessing, with re-customized options (here: 'events', to define how trials are to be
% extracted). We assume that the data contains events of the type 'user-action' whenever the user
% performs a (not further known) action, such as one of the imagined hand gestures.
%  
%   testdata = io_loadset('data sets/john/unknown_gestures.eeg')
%   testdataproc = bci_preproc(testdata, model, 'events',{'user-action'})
%   predictions = bci_predict(model, testdataproc, 'process',0);
%
%   disp('the user performed the following actions: ' num2str(predictions{2}(argmax(predictions{3}')))]);
%
% In this case, we obtain predictions for every trial of the data (this time relative to events of
% the type 'user-action'), whereas the the predictions will still be either 1,2, or 3, since the
% predictive part of the model has not been changed. Further, since the activity of the user was not
% known a priori, a meaningful loss cannot be computed. Note that if the preprocessing of the model
% contains filters that are adaptively tuned to the data (such as ICA), this strategy will give
% unexpected results, because bci_preproc re-runs the preprocessing (i.e. readapts the filters) for
% the test data, which gives intermediate values that are not expected by the rest of the model.
%
% For these reasons, the recommended way to obtain predictions for a data set at the times of
% arbitrary events, or in arbitrary intervals, is to use onl_stream, which streams the given raw
% data set into the model and obtains predictions at the desired time points. bci_predict should
% only be used to evaluate the performance of models on data set where the desired outcomes are
% known and encoded in the same way as in the respective calibration data set of the used model
% (i.e. for pure offline analyses). See onl_stream for further details.
%
%
% In:
%   Model  : a predictive model, as produced by bci_train
%
%   Data   : a raw data set from which to derive predictions; may also be a stream bundle
%
%   TargetMarkers : Target markers. By default as in the Model. Otherwise, this is a list of types
%                   of those markers around which data shall be used for BCI calibration; each
%                   marker type encodes a different target class (i.e. desired output value) to be
%                   learned by the resulting BCI model. See also help of bci_train.
%
%   EvaluationMetric : Evaluation metric. The metric to use in the assessment of model performance
%                      (via cross-validation). Can be empty, a string, or a function handle.
%                      See ml_calcloss() for the options (default: [] = auto-select between 
%                      kullback-leibler divergence ('kld'), mean square error ('mse'), mis-classification 
%                      rate ('mcr') and negative log-likelihood ('nll') depending on the type of the 
%                      target and prediction variables, further detailed in ml_calcloss())
%
%   Format       : format of the prediction; see utl_formatprediction (default: 'raw')
%
%   EventField : Event field to search for target markers. By default as in the Model.
%
%
% Out:
%   Prediction : the target variable as predicted by the model for every index in the data set
%
%   Loss       : a measure of the average loss of the model, w.r.t. to the target variable (as
%                evaluation; the allowed format is anything)
%
%   Stats      : additional statistics, as produced by the metric
%
%   Target     : original target variable, as determined by the target function
%
% Examples:
%   % given a predictive model and a continuous data set with markers that are compatible to those 
%   % that were used for training, derive per-epoch/per-marker BCI outputsm and estimate the loss 
%   % plus additional statistics
%   [predictions,loss,stats] = bci_predict(model,data)
%
%   % as before, but use a custom metric (here: mean-square error)
%   [predictions,loss,stats] = bci_predict('Data',data, 'Model',model, 'EvaluationMetric','mse')
%
% See also:
%   bci_predict, bci_annotate, onl_newpredictor, io_loadset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-24
dp;

% read arguments
opts = arg_define([0 2],varargin, ...
    arg_norep({'model','Model'},mandatory,[],'Predictive model. This is a model as previously computed via bci_train.'), ...
    arg_norep({'dataset','Data','data'},mandatory,[],'Data set. EEGLAB data set or stream bundle to use for prediction.'), ...
    arg({'markers','TargetMarkers'},'frommodel',[],'Assumed target markers. List of types of those markers around which data shall be used for BCI calibration; each marker type encodes a different target class (i.e. desired output value) to be learned by the resulting BCI model. This can be specified either as a cell array of marker-value pairs, in which case each marker type of BCI interest is associated with a particular BCI output value (e.g., -1/+1), or as a cell array of marker types (in which case each marker will be associated with its respective index as corresponding BCI output value, while nested cell arrays are also allowed to group markers that correspond to the same output value). See help of set_targetmarkers for further explanation.','type','expression'), ...
    arg({'metric','EvaluationMetric','Metric'},'auto',{'auto','mcr','mse','smse','nll','kld','mae','max','rms','bias','medse','auc','cond_entropy','cross_entropy','f_measure'},'Evaluation metric. The metric to use in the assessment of model performance (via cross-validation); see also ml_calcloss.'), ...
    arg({'outformat','Format','format'},'raw',{'raw','expectation','distribution','mode'},'Prediction format. See utl_formatprediction.'),...
    arg({'field','EventField'}, 'frommodel', [], 'Assumed target field. Event field to search for target markers, provided as a string.','type','char','shape','row'));

% input validation
if ~isstruct(opts.model) || ~isscalar(opts.model)
    error('The given Model argument must be a 1x1 struct.'); end
if ~isfield(opts.model,'tracking') || ~all(isfield(opts.model.tracking,{'prediction_function','filter_graph','prediction_channels'}))
    error('The given Model argument is lacking some required fields (required are: .tracking.prediction_function, .tracking.filter_graph, .tracking.prediction_channels), but got: %s',hlp_tostring(opts.model,10000)); end

% evaluate and check the prediction funtion
if ischar(opts.model.tracking.prediction_function)
    % prediction function given as a string
    if strncmp(opts.model.tracking.prediction_function,'Paradigm',8)
        % class reference: instantiate
        try
            instance = eval(opts.model.tracking.prediction_function); %#ok<NASGU>
        catch e
            error('Failed to instantiate paradigm class (%s) with error: %s',opts.model.tracking.prediction_function,e.message);
        end
        opts.model.tracking.prediction_function = eval('@instance.predict');
    else
        % some other function
        if ~exist(opts.model.tracking.prediction_function,'file')
            error('The given prediction function was not found on the MATLAB path: %s',opts.model.tracking.prediction_function); end
        opts.model.tracking.prediction_function = str2func(opts.model.tracking.prediction_function);
    end
end

% use the model's target field if desired
if strcmp(opts.field,'frommodel')
    opts.field = opts.model.control_options.field; end

% use the model's target markers if desired
if strcmp(opts.markers,'frommodel')
    opts.markers = opts.model.control_options.markers; end

% string arguments are considered to be variants of the default metric
if isempty(opts.metric) || ischar(opts.metric) || (iscell(opts.metric) && all(cellfun(@ischar,opts.metric)))
    opts.metric = @(T,P)ml_calcloss(opts.metric,T,P); end
if ~has_stats(opts.metric)
    opts.metric = @(T,P)add_stats(opts.metric(T,P)); end

% uniformize data
dataset = opts.dataset;
if iscell(dataset)
    error('The bci_predict function cannot be applied to dataset collections -- you need to apply it to each dataset individually.'); end
if ~isfield(dataset,'streams')
    dataset = struct('streams',{{dataset}}); end
% and annotate with target markers (there are 3 possible TargetMarker formats to handle)
if length(opts.markers) == 1 && ischar(opts.markers{1}) && strcmp(opts.markers{1}, 'actualvalues')
    dataset.streams{1} = set_targetmarkers('Signal',dataset.streams{1},'EventMap',opts.markers,'EpochBounds',opts.model.epoch_bounds, 'EventField', opts.field );
elseif all(cellfun('isclass',opts.markers,'char') | cellfun('isclass',opts.markers,'cell'))
    dataset.streams{1} = set_targetmarkers('Signal',dataset.streams{1},'EventTypes',opts.markers,'EpochBounds',opts.model.epoch_bounds, 'EventField', opts.field);
else
   dataset.streams{1} = set_targetmarkers('Signal',dataset.streams{1},'EventMap',opts.markers,'EpochBounds',opts.model.epoch_bounds, 'EventField', opts.field );
end
% and do final checks and fixups
dataset = utl_check_bundle(dataset);

% attach the passed stream bundle to the model's filter graph
model = opts.model;
resolved_graph = utl_resolve_streams(model.tracking.filter_graph,dataset,model.tracking.prediction_channels);
% process each chain of the filter graph
dataset.streams = cell(1,length(resolved_graph));
for g=1:length(resolved_graph)
    dataset.streams{g} = exp_eval_optimized(resolved_graph{g}); end
% apply the prediction function to the data
prediction = model.tracking.prediction_function(dataset,model);

try
    % try to derive the target variable
    target = set_gettarget(dataset);
catch e
    % a target variable is not necessarily present
    disp('NOTE: Did not find a target variable in the given data (' + e.message + ')');
    target = [];
end
if ~isempty(target)
    % and evaluate the performance, given the metric
    [measure,stats] = opts.metric(target,prediction);
else
    measure = NaN;
    stats = struct();
end
stats.target = target;
stats.prediction = prediction;
stats.is_result = true;
stats.timestamp = now;

% finally format the predictions
prediction = utl_formatprediction(prediction,opts.outformat);


% test whether the given metric supplies stats or not
function result = has_stats(metric)
try 
    [x,y] = metric([],[]);  %#ok<NASGU,ASGLU>
    result = true;
catch e
    result = ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs','MATLAB:unassignedOutputs'})); 
end

% add stats to the result of some metric
function [result,stats] = add_stats(result)
stats.value = result;
