function [measure,stats] = utl_nested_crossval(data, varargin)
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

% get options struct & assign some defaults as documented
opts = hlp_varargin2struct(varargin, ...
    ... % arguments that specifically control the outer cross-validation
    'eval_scheme',10, ...
    ... % arguments that control the grid search
    'opt_scheme',5, ...
    'argform','clauses', ...
    'args',{});

% if parameters are to be searched
if is_needing_search(opts.argform,opts.args)
    % use the utl_searchmodel() meta-function as the trainer in the subsequent regular CV
    % with opts.opt_scheme as the (nested) search scheme
    opts.trainer = @(P,varargin) utl_searchmodel(P,opts,'scheme',opts.opt_scheme); 
end

% run a cross-validation, with opts.eval_scheme as the search scheme
[measure,stats] = utl_crossval(data, opts, 'scheme',opts.eval_scheme);
