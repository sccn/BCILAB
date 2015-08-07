function [measure,stats] = utl_nested_crossval(varargin)
% Run a generic nested cross-validation over indexable data.
% [Measure, Stats] = utl_nested_crossval(Data, Arguments...)
%
% In:
%   Data        :   some data that can be partitioned using index sets
%
%   Arguments   :   mandatory arguments:
%                   'trainer': training function; receives a partition of the data (as produced by the partitioner),
%                              possibly some further arguments as specified in args, and returns a model (of any kind)
%
%                   'tester' : testing function; receives a model (as produced by the trainer) and a partition of the data
%                              (as produced by the partitioner), and returns a prediction for every index of the input data,
%                              in a format allowed by ml_predict
%
%                   'args': argument (ranges) for the training function, cell array (may be empty)
%                           specified as in utl_gridsearch (format controlled via argform)
%
%                   optional arguments:
%                   'opt_scheme': cross-validation scheme used while searching for the optimal model (default: 5)
%                                  uses the scheme format of utl_crossval(), see below
%                   'eval_scheme': cross-validation scheme used while evaluating the optimal model's performance (default: 10)
%                                  uses the scheme format of utl_crossval():
%                                   * 0: skip CV, return NaN and empty statistics
%                                   * k: k-fold randomized CV (or, if 0<k<1, k-holdhout CV)
%                                   * [r k]: r-time repartitioned k-fold randomized CV (or, if 0<k<1, k-holdout CV (with k a fraction))
%                                   * 'loo': leave-one-out CV
%                                   * {'chron', k} or {'block', k}: k-fold chronological/blockwise CV
%                                   * {'chron', k, m} or {'block', k, m}: k-fold chronological/blockwise CV with m indices margin width (between training and test set)
%
%                   'partitioner': partitioning function for the data, receives three parameters: (data, index vector, packed trainer args OR model)
%                                   * if the index vector is empty, should return the highest index in the data
%                                   * otherwise, it should return data subindexed by the index vector
%                                  default: provides support for cell arrays, numeric arrays, struct arrays and {Data,Target} cell arrays
%                                  note: the third parameter is for convenience and may optionally be taken into account; trainer args are passed packed into 
%                                  a cell array for both index set generation and computation of the training partition(s), and the model is passed for the
%                                  computation of testing partitions
%
%                   'target': a function to derive the target variable from a partition of the data (as produced by the partitioner),
%                             for evaluation; the allowed format is anything that may be output by ml_predict
%                             default: provides support for {Data,Target} cell arrays
%
%                   'metric': metric to be employed, applied both in each fold and the aggregated data over all folds
%                              * function handle: a custom, user-supplied loss function; receives target data in the first 
%                                argument and prediction data in the second argument; each can be in any format that can be 
%                                produced by ml_predict (but can be expected to be mutually consistent).
%                                shall return a real number indicating the summary metric over all data, and optionally 
%                                additional statistics in a struct
%                              * string: use ml_calcloss, with 'metric' determining the loss type
%                              * default: use 'mcr','mse','nll','kld', depending on supplied target & prediction data formats
%
%                   'argform': format of the argument ranges, either 'direct' or 'clauses' (default: 'clauses')
%                               'direct': search ranges are directly specified as arrays
%                               'clauses': search ranges are specified using the search() clause
%                   
%                   further arguments (same as in utl_crossval)
%                   'repeatable': whether the randomization procedure shall give repeatable results 
%                                 (default: 1); different numbers (aside from 0) give different
%                                 repeatable runs, i.e. the value determines the randseed
%
%                   'forcestats': enforce a nested cross-validation in order to obtain stats, even
%                                 when no parameter search is necessary (default: 0)
%
%                   'collect_models': whether to return models trained for each fold (default: false)
%
%                   'cache_fold_results' : whether to cache the per-fold results. (default: false)
%
%                   'only_cached_results' : load only results that are in the cache. (default: false)
%
%                   'tolerate_exceptions' : tolerate exceptions during training step. (default: false)
%
%                   further arguments (same as in par_beginschedule, listed below for convenience)
%                   'engine_ps': the parallelization engine to be used for the parameter search (default: 'local')
%
%                   'engine_cv': the parallelization engine to be used for the cross-validation (default: 'global')
%
%                   'pool'     : node pool to be used for parallelization, when using the BLS scheduler (default: 'global')
%
%                   'policy'   : scheduling policy to be used, when using the BLS scheduler (default: 'global')
%
% Out:
%   Measure : a measure of the overall performance of the trainer/tester combination, w.r.t. to the target variable returned by gettarget,
%             computed by the metric
%
%   Stats   : additional statistics, as produced by the metric
%
% Example:
%   % assuming a feature matrix called trials and a label vector called targets, sized as:
%   %  trials: [NxF] array of training instances
%   %  targets: [Nx1] array of labels
%
%   % cross-validate using sparse logistic regression
%   [loss,stats] = utl_nested_crossval({trials,targets}, 'args',{{'logreg' search(2.^(-6:1:10)) 'variant' 'l1'}})  
%
% See also:
%   utl_crossval, utl_searchmodel
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22
dp;

% parse options
opts = arg_define(0:1,varargin, ...
    ... % training data
    arg_norep({'data','Data'},[],[],'Data that can be partitioned using index sets. Such as, for example, {X,y} with X being the [NxF] training data and y being the [Nx1] labels (N=#trials, F=#features).'), ...
    ... % method to use and its arguments
    arg({'args','TrainingArguments'},{},[],'Arguments to the training function. Packed in a cell array. Search ranges are defined according to ArgumentFormat; see also utl_gridsearch for additional clarification on this. Note: if using the default trainer/tester functions, args must at least specify the learning function to be used, optionally followed by arguments to that learning function, e.g. {''lda''} for linear discriminant analysis or {''logreg'',''lambda'',0.1} for logistic regression with ''lambda'',0.1 passed in as user parameters (see ml_train* functions for options).','type','expression'), ...
    arg({'trainer','TrainingFunction'},'@ml_train',[],'Training function. Receives a partition of the data (as produced by the specified partitioner), possibly some further arguments as specified in args, and returns a model (of any kind).','type','expression'), ...
    arg({'tester','PredictionFunction'},'@utl_default_predict',[],'Prediction function. Receives a partition of the data (as produced by the partitioner) and a model (as produced by the trainer), and returns a prediction for every index of the input data, in one of the output formats permitted by ml_predict.','type','expression'), ...
    ... % nested cross-validation parameters
    arg({'eval_scheme','EvaluationScheme'},10,[],'Outer cross-validation scheme. Defines what parts of the data shall be used for training or testing in each fold. Supports many different formats, see documentation.','type','expression'), ...
    arg({'opt_scheme','OptimizationScheme'},5,[],'Parameter optimization CV scheme. Defines what parts of the data shall be used for training or testing in each fold. Supports many different formats, see documentation.','type','expression'), ...
    arg({'metric','EvaluationMetric'},'auto',{'auto','mcr','auc','mse','medse','smse','mae','medae','smae','smedae','max','sign','rms','bias','kld','nll','cond_entropy','cross_entropy','f_measure'},'Evaluation metric. Loss measure employed to measure the quality of predictions on a test partition given its known target values; applied both to results in each fold and results aggregated over all folds. Can be either a predefined metric supported by ml_calcloss, or a custom function handle which receives an array of target values in the first argument and an array of predictions in the second argument; each can be in any format that can be produced by ml_predict (but can be expected to be mutually consistent). Shall return a real number indicating the summary metric over all data, and optionally additional statistics in a second output struct.','typecheck',false), ...
    arg({'repeatable','RepeatableResults'},1,[],'Produce repeatable (vs. randomized) results. If nonzero, the value is taken as the random seed to use for the performed calculation. This way, different numbers (aside from 0) give different repeatable runs.'), ...
    ... % support for custom data representations
    arg({'target','TargetFunction'},'@utl_default_target',[],'Target extraction function. A function to derive the target variable from a partition of the data (as produced by the partitioner), for evaluation; the allowed format is anything that may be output by ml_predict.','type','expression'), ...
    arg({'partitioner','PartitioningFunction'},'@utl_default_partitioner',[],'Partitioning function. See documentation for the function contract.','type','expression'), ...
    ... % parallel computing settings
    arg({'engine_gs','ParallelEngineGS'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine for grid search. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'engine_cv','ParallelEngineCV'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine for outer cross-validation. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'engine_ncv','ParallelEngineNCV'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine for inner (nested) cross-validation. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'pool','WorkerPool'},'global',[], 'Worker pool to use. This is typically a cell array, but can also be the string ''gobal'', which stands for the currently globally set up worker pool (see global tracking variable).','type','expression'), ...
    arg({'policy','ReschedulingPolicy'},'global',[], 'Rescheduling policy. This is the name of the rescheduling policy function that controls if and when tasks are being rescheduled. If set to global, the current global setting will be used.'), ...    
    ... % misc arguments    
    arg({'argform','ArgumentFormat'},'clauses',{'clauses','direct'},'Argument search range format. If set to clauses, search ranges for individual arguments can be given using the search() expression; if set to direct, a cell array is expected each of whose elements is a cell array that represents a particular parameter set to try for the training function.'), ...
    arg({'forcestats','ForceStatistics'},0,[],'Enforce a cross-validation to obtain stats. Even when no search is required.'), ...
    arg({'cache_fold_results','CacheFoldResults'},false,[],'Whether to cache the per-fold results. This is meant to be used when running very long-running computations on machines that crash frequently enough that partial results need to be saved. In this case, any previously computed results will be loaded from disk.'), ...
    arg({'only_cached_results','OnlyCachedResults'},false,[],'Load only results that are in the cache. This will not run any computations (aside from pre-checks, that can be disabled by setting NoPrechecks to true).'), ...
    arg({'no_prechecks','NoPrechecks'},false,[],'Skip pre-checks that access the data. This can save some time when it would take very long to load the data, especially when performing parallel computation.'), ...
    arg({'tolerate_exceptions','TolerateExceptions'},false,[],'Tolerate exceptions during training. If this happens, folds where the training function yielded errors will be skipped.'), ...
    arg({'collect_models','CollectModels'},false,[],'Collect models per fold. Note that this increases the amount of data returned.'));

% if parameters are to be searched...
if is_needing_search(opts.argform,opts.args)
    % use the utl_searchmodel() meta-function as the trainer in the subsequent regular CV
    % with opts.opt_scheme as the (nested) search scheme
    opts.trainer = @(data,varargin) utl_searchmodel(data,rmfield(opts,{'eval_scheme','opt_scheme','engine_cv','data'}),'scheme',opts.opt_scheme); 
end

% run a cross-validation, with opts.eval_scheme as the search scheme
[measure,stats] = utl_crossval(rmfield(opts,{'eval_scheme','opt_scheme','engine_gs','engine_ncv','argform','forcestats'}), 'scheme',opts.eval_scheme);
