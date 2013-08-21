function [model,stats] = utl_searchmodel(data, varargin)
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

% get options struct, and impose documented defaults (any further arguments
% have no defaults)
opts = hlp_varargin2struct(varargin, ...
    ... % arguments for the model search process
    'trainer',@ml_train, ...    % the default training function
    'forcestats',0, ...         % whether a cross-validation should always be done
    ... % arguments that control the inner cross-validation
    'engine_ncv','global', ...  % the engine_ncv argument will be used as the engine_cv argument of utl_crossval
    ... % arguents that control the grid search
    'engine_gs','global', ...   % the grid searching engine
    'argform','clauses', ...    % by default we use the 'clauses' grid search format
    'args',{});
    
% if there are parameters to be searched...
if opts.forcestats || is_needing_search(opts.argform,opts.args)
    
    % set up the appropriate options struct that will be used for the nested 
    % cross-validation; in particular, all our options are forwarded to it, 
    % but the engine_gs parameter is peeled off...
    nestedcv_ctrl = rmfield(opts,'engine_gs');
    % ... and the value for engine_cv (whether present or not) is substituted by the value of engine_ncv
    nestedcv_ctrl.engine_cv = opts.engine_ncv;

    % set up the objective function that will be used for the grid search,
    % which is utl_crossval, using the control options defined above
    % The objective function's free parameters are passed down into
    % the cross-validation as its 'args' parameter, which is the arguments
    % that are forwarded by it to the 'trainer' function
    objective_function = @(varargin) utl_crossval(data,nestedcv_ctrl,'args',varargin);
    
    % set up the appropriate options struct for utl_gridsearch; in particular, 
    % the parallelization & clauses behavior of or the gridsearch is determined according to the opts...
    gridsearch_ctrl = hlp_struct2varargin(opts,'restrict',{'argform','pool','policy','engine_gs'});
    % ... and the objective function ('func') is the one defined above 
    gridsearch_ctrl = [gridsearch_ctrl {'func',objective_function}];
    
    % now execute the grid search over every specified combination in opts.args
    [stats.bestidx,stats.inputs,stats.outputs] = utl_gridsearch(gridsearch_ctrl, opts.args{:});
    
else
    
    % otherwise produce dummy stats if neither necessary nor forced
    stats = struct('bestidx',1, 'inputs',{{opts.args}}, 'outputs',{{NaN,struct()}});
    
end

% using the best args, send the complete data through the training function & compute a model on it
bestargs = stats.inputs{stats.bestidx};
model = opts.trainer(data,bestargs{:});
