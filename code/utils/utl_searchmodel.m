function [model,stats] = utl_searchmodel(varargin)
% Find the best predictive model out of a parameterized set, via cross-validation.
% [Model,Stats] = utl_searchmodel(Data, Arguments...)
%
% When a predictive model is learned from some data, several of its key parameters can usually be 
% efficiently optimized given the data (namely, if their solutions are available in closed form, or 
% if they can be found by solving a well-behaved optimization problem). Some parameters, however,
% usually remain, which can not be optimized by the given learning function. If their optimal values
% are unknown but assumed to lie in a small set of possibilities, a general-purpose approach to find 
% them as well, is exhaustive search over all possibilities, coupled with a cross-validation to 
% assess the predictive performance for each possibility on the available data. 
%
% This function performs this task, and can jointly optimize all unknown parameters of a learning
% function. It requires a data set and (because a cross-validation is implicitly performed) all
% parameters that would be required to apply utl_crossval to the given data. The function whose
% parameters shall be optimized is the 'trainer' (i.e., model learning) function. For any of its
% parameters which has multiple possibilities (over which to optimize), these possibilities can be
% specified using the same syntax as in the grid searching function utl_gridsearch (where the 
% default here is to specify them via search() clauses).
%
% In:
%   Data :   some data that can be partitioned using index sets
%
%   Arguments : mandatory arguments:
%               'trainer': training function; receives a partition of the data (as produced by 
%                          the partitioner), some further arguments (as specified in args), and
%                          returns a model (of any kind)
%
%               'tester' : testing function; receives a model (as produced by the trainer) and a 
%                          partition of the data (as produced by the partitioner), and returns a
%                          prediction for every index of the input data, in a format allowed by
%                          ml_predict
%
%               'args': arguments to the training function, cell array (default: empty)
%                       specified as in utl_gridsearch (format controlled via argform)
%
%               optional arguments (same as in utl_crossval, listed below for convenience)
%               'scheme': cross-validation scheme, can be one of the following: (default: 10)
%                 * 0: skip CV, return NaN and empty statistics
%                 * k: k-fold randomized CV (or, if 0<k<1, k-holdhout CV)
%                 * [r k]: r-time repartitioned k-fold randomized CV (or, if 0<k<1, 
%                   k-holdout CV (with k a fraction))
%                 * 'loo': leave-one-out CV
%                 * {'chron', k} or {'block', k}: k-fold chronological/blockwise CV
%                 * {'chron', k, m} or {'block', k, m}: k-fold chronological/blockwise CV with 
%                   m indices margin width (between training and test set)
%
%               'partitioner': partitioning function for the data, receives three parameters: 
%                              (data, index vector, packed trainer args OR model)
%                               * if the index vector is empty, should return the highest index 
%                                 in the data
%                               * otherwise, it should return data subindexed by the index vector
%                              default: provides support for cell arrays, numeric arrays, struct 
%                                       arrays and {Data,Target} cell arrays
%                              note: the third parameter is for convenience and may optionally be 
%                                    taken into account; trainer args are passed packed into a cell
%                                    array for both index set generation and computation of the
%                                    training partition(s), and the model is passed for the
%                                    computation of testing partitions
%
%               'target': a function to derive the target variable from a partition of the data (as 
%                         produced by the partitioner), for evaluation; the allowed format is
%                         anything that may be output by ml_predict default: provides support for
%                         {Data,Target} cell arrays
%
%               'metric': metric to be employed, applied both in each fold and the aggregated data 
%                         over all folds
%                          * function handle: a custom, user-supplied loss function; receives target 
%                            data in the first argument and prediction data in the second argument;
%                            each can be in any format that can be produced by ml_predict (but can
%                            be expected to be mutually consistent). shall return a real number
%                            indicating the summary metric over all data, and optionally additional
%                            statistics in a struct
%                          * string: use ml_calcloss, with 'metric' determining the loss type
%                          * default: use 'mcr','mse','nll','kld', depending on supplied target & 
%                            prediction data formats
%
%               'argform': format of the argument ranges, either 'direct' or 'clauses'
%                          (default: 'clauses')
%                           'direct': search ranges are directly specified as arrays
%                           'clauses': search ranges are specified using the search() clause
%
%               'forcestats': enforce a cross-validation to obtain stats, even when no search is 
%                             necessary (default: 0)
%
%               further arguments (same as in utl_crossval)
%               'cache_fold_results' : whether to cache the per-fold results. (default: false)
%
%               'only_cached_results' : load only results that are in the cache. (default: false)
%
%               'tolerate_exceptions' : tolerate exceptions during training step. (default: false)
%
%               further arguments (same as in par_beginschedule, listed below for convenience)
%               'engine_gs': the parallelization engine to be used for the grid search
%                            (default: 'local')
%
%               'engine_ncv': the parallelization engine to be used for the nested cross-validation
%                             (default: 'global')
%
%               'pool'     : node pool to be used for parallelization, when using the BLS scheduler
%                            (default: 'global')
%
%               'policy'   : scheduling policy to be used, when using the BLS scheduler
%                            (default: 'global')
%
% Out:
%   Model : the best model, as determined through cross-validation and and grid search 
%           w.r.t. to a measure of the overall performance of the trainer/tester combination
%
%   Stats : additional statistics over all tested argument combinations, as produced by the metric
%
% See also:
%   utl_crossval, utl_gridsearch, utl_nested_crossval
%
% Example:
%   % assuming a feature matrix called trials and a label vector called targets, sized as:
%   %  trials: [NxF] array of training instances
%   %  targets: [Nx1] array of labels
%
%   % find the best-performing SVM classifier on the data (here: ignoring the kernel scale parameter)
%   [model,stats] = utl_searchmodel({trials,targets}, 'args',{{'svm' search(2.^(-5:2:15))}}) 
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-22
dp;

% parse arguments
opts = arg_define(0:1,varargin, ...
    ... % training data
    arg_norep({'data','Data'},[],[],'Data that can be partitioned using index sets. Such as, for example, {X,y} with X being the [NxF] training data and y being the [Nx1] labels (N=#trials, F=#features).'), ...
    ... % method to use and its arguments
    arg({'args','TrainingArguments'},{},[],'Arguments to the training function. Packed in a cell array. Search ranges are defined according to ArgumentFormat; see also utl_gridsearch for additional clarification on this. Note: if using the default trainer/tester functions, args must at least specify the learning function to be used, optionally followed by arguments to that learning function, e.g. {''lda''} for linear discriminant analysis or {''logreg'',''lambda'',0.1} for logistic regression with ''lambda'',0.1 passed in as user parameters (see ml_train* functions for options).','type','expression'), ...
    arg({'trainer','TrainingFunction'},'@ml_train',[],'Training function. Receives a partition of the data (as produced by the specified partitioner), possibly some further arguments as specified in args, and returns a model (of any kind).','type','expression'), ...
    arg({'tester','PredictionFunction'},'@utl_default_predict',[],'Prediction function. Receives a partition of the data (as produced by the partitioner) and a model (as produced by the trainer), and returns a prediction for every index of the input data, in one of the output formats permitted by ml_predict.','type','expression'), ...
    ... % nested cross-validation parameters
    arg({'scheme','EvaluationScheme'},10,[],'Cross-validation scheme. Defines what parts of the data shall be used for training or testing in each fold. Supports many different formats, see documentation.','type','expression'), ...
    arg({'metric','EvaluationMetric'},'auto',{'auto','mcr','auc','mse','medse','smse','mae','medae','smae','smedae','max','sign','rms','bias','kld','nll','cond_entropy','cross_entropy','f_measure'},'Evaluation metric. Loss measure employed to measure the quality of predictions on a test partition given its known target values; applied both to results in each fold and results aggregated over all folds. Can be either a predefined metric supported by ml_calcloss, or a custom function handle which receives an array of target values in the first argument and an array of predictions in the second argument; each can be in any format that can be produced by ml_predict (but can be expected to be mutually consistent). Shall return a real number indicating the summary metric over all data, and optionally additional statistics in a second output struct.','typecheck',false), ...
    arg({'repeatable','RepeatableResults'},1,[],'Produce repeatable (vs. randomized) results. If nonzero, the value is taken as the random seed to use for the performed calculation. This way, different numbers (aside from 0) give different repeatable runs.'), ...
    ... % support for custom data representations
    arg({'target','TargetFunction'},'@utl_default_target',[],'Target extraction function. A function to derive the target variable from a partition of the data (as produced by the partitioner), for evaluation; the allowed format is anything that may be output by ml_predict.','type','expression'), ...
    arg({'partitioner','PartitioningFunction'},'@utl_default_partitioner',[],'Partitioning function. See documentation for the function contract.','type','expression'), ...
    ... % parallel computing settings
    arg({'engine_gs','ParallelEngineGS'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine for grid search. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'engine_ncv','ParallelEngineNCV'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine for cross-validation. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
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

% if there are parameters to be searched...
if opts.forcestats || is_needing_search(opts.argform,opts.args)
    
    % the nested cross-validation (utl_crossval) receives a specific subset of our arguments
    nestedcv_opts = hlp_struct2varargin(opts,'suppress',{'engine_gs','forcestats','args','argform'},'rewrite',{'engine_ncv','engine_cv'});

    % the objective function of our search is utl_crossval with the above-defined options; plus, it
    % has free arguments that are passed into the cross-validation as the 'args' parameter (i.e.,
    % the arguments that it passes into its 'trainer' function)
    objfun = @(varargin) utl_crossval(nestedcv_opts{:}, 'args',varargin);
    
    % the grid search receives the above-defined objective function and a subset of our arguments
    search_opts = {'argform',opts.argform, 'engine_gs',opts.engine_gs, 'policy',opts.policy, ...
        'pool',opts.pool, 'func',objfun};
    
    % now execute the grid search over every argument combination specified in opts.args
    [stats.bestidx,stats.inputs,stats.outputs] = utl_gridsearch(search_opts, opts.args{:});
    
else
    
    % otherwise produce dummy stats if neither necessary nor forced
    stats = struct('bestidx',1, 'inputs',{{opts.args}}, 'outputs',{{NaN,struct()}});
    
end

% finally, using the best args, send the input data through the training function & compute a model
bestargs = stats.inputs{stats.bestidx};
model = opts.trainer(opts.data,bestargs{:});
