function [measure,stats] = utl_crossval(varargin)
% Run a generic cross-validation over indexable data.
% [Measure, Stats] = utl_crossval(Data, Arguments...)
%
% Cross-validation [1] is a data resampling technique in which per each iteration (or "fold"), a
% model is formed given a subset of the data (called the training set), and then its quality is
% tested on the remaining portion of the data (called the test set). Most applications of
% cross-validation ensure that each observation (or trial) in the data is used for testing in at
% least one fold. The most common version is k-fold cross-validation, where the data is split into k
% parts, and throughout k iterations of the algorithm, each of the parts is chosen as the test set
% exactly once, while the remaining k-1 parts are used jointly to train the model that shall be
% tested. Thus, each part of the data is used in k-1 training sets.
%
% The partitions of the data are traditionally chosen in a randomized fashion, under the assumption
% that each trial is identically and independently distributed w.r.t. the others. This assumption
% is wrong in time-series data (i.e. in almost all BCI cases). For this reason, the recommended 
% cross-validation scheme in BCILAB is "chronological" (a.k.a. "block-wise") cross-validation, where
% the trials within the test set are taken from a compact interval of the underyling time series.
% 
% In addition, it is recommended to exclude trials that are close to any of the test set trials 
% from use in the respective training set (these excluded trials are called safety margins).
%
% The utl_crossval function is primarily used for testing predictive models, i.e., models which are
% learned to be able to estimate some "target" value for each data trial. Thus, the quality of a
% model given some test data is primarily assessed by letting it predict target values for each test
% trial, and comparing the outcomes with the known values for the given set of test trials.
%
% The data to be cross-validated can be in a variety of formats (as long as utl_crossval can
% determine how to partition it), and a custom partitioning function can be supplied for completely
% free-form data (such as EEGLAB data sets, which are handled by utl_dataset_partitioner). Aside
% from the data, usually two functions need to be passed - one to learn a model from some data
% portion (the 'trainer'), and one to test the learned model on some data portion (the 'tester'), by
% making predictions. Further customization options include choice of the comparison metric (or 
% 'loss' function), choice of the function which extracts known target values from the data (as the 
% ground truth used in the comparison), as well as cluster parallelization options.
%
% In:
%   Data : some data that can be partitioned using index sets, such as, for example, {X,y} with X 
%          being the [NxF] training data and y being the [Nx1] labels (N=#trials, F=#features). Any
%          data format that is supported by the given trainer, tester, partitioner and target
%          functions is supported by utl_crossval.
%
%   Arguments : optional name-value pairs specifying the arguments:
%               'trainer': training function; receives a partition of the data (as produced by 
%                          the specified partitioner), possibly some further arguments as
%                          specified in args, and returns a model (of any kind) (default:
%                          @ml_train; natively supports the {X,y} data format)
%
%               'tester' : testing function; receives a partition of the data (as produced by 
%                          the partitioner) and a model (as produced by the trainer), and
%                          returns a prediction for every index of the input data, in one of the
%                          output formats permitted by ml_predict (default: @ml_predict)
%
%               'scheme': cross-validation scheme, can be one of the following formats (default: 10)
%                Generally, in the following k is the number of folds in k-fold CV, or if 0<k<1 
%                it is the fraction of test trials per fold, as in p-holdout CV. Also m is the
%                number of trials of margin between any training and test trial, or if 0<m<1, it is
%                given as a fraction of the number of total trials in the data. r is the number of
%                repeats for the CV procedure (as in r times repeated k-fold CV).
%                * 0: skip CV, return NaN and empty statistics
%                * k: k-fold randomized CV (or, if 0<k<1, randomized k-holdhout CV)
%                * [r k]: r times repartitioned k-fold randomized CV (or, if 0<k<1, randomized 
%                         k-holdout CV (with k a fraction))
%                * [r k m]: r times repartitioned k-fold randomized CV (or, if 0<k<1, k-holdout CV 
%                           (with k a fraction)) with m indices margin width (between training and
%                           test set)
%                * 'loo': leave-one-out CV
%                * {'chron', k} or {'block', k}: k-fold chronological/blockwise CV (or k-holdout 
%                                                where the held-out test trials are at the end of
%                                                the data, if 0<k<1)
%                * {'chron', k, m} or {'block', k, m}: k-fold chronological/blockwise CV with m 
%                                                      indices margin width (between training and 
%                                                      test set) or k-holdout as above
%                * 'trainerr': Return the training-set error (a measure of the separability of the 
%                              data, not of the generalization ability of the classifier); it is 
%                              usually an error to report this number in a paper.
%
%               'partitioner': partitioning function for the data, receives two parameters: 
%                              (data, index vector)
%                              * if the index vector is empty: in the simplest case the function
%                                should return the highest index in the data. For more custom
%                                control it may return a cell array of the form {{train-partition-1,
%                                test-partition-1}, {train-partition-2, test-partition-2}, ...}. In
%                                this case it fully controls the cross-validation scheme for each
%                                fold; the *-partition-n arrays should ideally be index ranges, but
%                                they will not be inspected by utl_crossval; instead, to later
%                                partition the data for each fold, the respective *-partition-n
%                                array will be passed back in to the partitioner as the "index
%                                vector" argument (besides with the whole data), which allows the
%                                partitioner to decide its own format. In the most advanced case the
%                                per-fold cell arrays returned in response to a call with index
%                                vector may have a third element as in {train-partition-1,
%                                test-partition-1, train-options}; then it will be taken as a cell
%                                array of options that are to be passed into the trainer on the
%                                respective fold as extra arguments (this allows the trainer to be
%                                aware of the kind of partition that it received).
%                              * if the index vector is nonempty: the function should return data
%                                subindexed by the index vector
%
%                              default: a function that can partition cell arrays, numeric arrays, struct
%                                       arrays, {Data,Target} cell arrays, and EEGLAB dataset structs
%
%               'target': a function to derive the target variable from a partition of the data 
%                         (as produced by the partitioner), for evaluation; the allowed format is
%                         anything that may be output by ml_predict; default: provides support for
%                         {Data,Target} cell arrays
%
%               'metric': loss metric employed to measure the quality of predictions on a test 
%                         partition given its known target values; applied both to results in each 
%                         fold and results aggregated over all folds. can be one of the following:
%                         * function handle: a custom, user-supplied loss function; receives an 
%                           array of target values in the first argument and an array of predictions 
%                           in the second argument; each can be in any format that can be produced 
%                           by ml_predict (but can be expected to be mutually consistent). shall 
%                           return a real number indicating the summary metric over all data, and 
%                           optionally additional statistics in a second output struct
%                         * string: use ml_calcloss, with 'metric' determining the loss type to use
%                         * default/empty/'auto': use 'mcr','mse','nll','kld' depending on supplied 
%                           target and prediction data formats; see also ml_calcloss
%
%               'args': optional arguments to the training function, packed in a cell array 
%                       (default: empty)
%                       note: if using the default trainer/tester functions, args must at least 
%                             specify the learning function to be used, optionally followed by 
%                             arguments to that learning function, e.g. {'lda'} for linear
%                             discriminant analysis or {'logreg','lambda',0.1} for logistic regression 
%                             with 'lambda',0.1 passed in as user parameters (see ml_train* functions 
%                             for options)
%
%               'repeatable': whether the randomization procedure shall give repeatable results 
%                             (default: 1); different numbers (aside from 0) give different
%                             repeatable runs, i.e. the value determines the randseed
%
%               'collect_models': whether to return models trained for each fold (default: false)
%
%               'engine_cv': parallelization engine to use (default: 'global'); see par_beginsschedule
%                   
%               'pool'  : worker pool to use (default: 'global'); see par_beginsschedule
%
%               'policy': scheduling policy to use (default: 'global'); see par_beginschedule
%
%               'cache_fold_results' : whether to cache the per-fold results. (default: false)
%
%               'only_cached_results' : load only results that are in the cache. (default: false)
%
%               'no_prechecks' : skip pre-checks that access the data. (default: false)
%
%               'collect_models' : collect models per fold. (default: false)
%
% Out:
%   Measure : a measure of the overall performance of the trainer/tester combination, w.r.t. to the 
%             target variable returned by the target function. Computed according to the selected metric.
%
%   Stats   : additional statistics, as produced by the selected metric
%
% Example:
%   % assuming a feature matrix called trials and a label vector called targets, sized as:
%   %  trials: [NxF] array of training instances
%   %  targets: [Nx1] array of labels
%
%   % cross-validate using (shrinkage) linear discriminant analysis
%   [loss,stats] = utl_nested_crossval({trials,targets}, 'args',{'lda'})  
%
%   % cross-validate using hierarchical kernel learning with a specific kernel
%   [loss,stats] = utl_nested_crossval({trials,targets}, 'args',{{'hkl' 'kernel' 'hermite'}})  
%
% Configuration Examples:
%   A simple training function would be:
%     @ml_train, with args being {'lda'} -- for this case, the data X (size NxF) and labels y 
%                                           (size Nx1) should be supplied as {X,y} to utl_crossval
%
%   A simple prediction function would be:
%     @ml_predict
%
%   A simple partitioner for epoched EEG data sets would be:
%     function result = my_partitioner(data,indices,misc)
%     if isempty(indices)
%         result = data.trials;
%     else
%         result = exp_eval(set_selepos(data,indices));
%     end
%
%   A simple mean-square error loss metric would be:
%     my_metric = @(target,predicted) mean((target-predicted).^2);
%
%   A simple target extraction function for epoched EEG data sets would be (assuming that there is 
%   an epoch-associated target value):
%     my_target = @(data) [data.epoch.target];
%
% References:
%   [1] Richard O. Duda, Peter E. Hart, David G. Stork, "Pattern Classification" 
%       Wiley Interscience, 2000
%
% See also:
%   utl_evaluate_fold, bci_train, utl_searchmodel, utl_nested_crossval
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-07
dp;

% read arguments
opts = arg_define(0:1,varargin, ...
    ... % evaluation data
    arg_norep({'data','Data'},[],[],'Data that can be partitioned using index sets. Such as, for example, {X,y} with X being the [NxF] training data and y being the [Nx1] labels (N=#trials, F=#features).'), ...
    ... % evaluation parameters
    arg({'scheme','EvaluationScheme'},10,[],'Cross-validation scheme. Defines what parts of the data shall be used for training or testing in each fold. Supports many different formats, see documentation.','type','expression'), ...
    arg({'metric','EvaluationMetric'},'auto',{'auto','mcr','auc','mse','medse','smse','mae','medae','smae','smedae','max','sign','rms','bias','kld','nll','cond_entropy','cross_entropy','f_measure'},'Evaluation metric. Loss measure employed to measure the quality of predictions on a test partition given its known target values; applied both to results in each fold and results aggregated over all folds. Can be either a predefined metric supported by ml_calcloss, or a custom function handle which receives an array of target values in the first argument and an array of predictions in the second argument; each can be in any format that can be produced by ml_predict (but can be expected to be mutually consistent). Shall return a real number indicating the summary metric over all data, and optionally additional statistics in a second output struct.','typecheck',false), ...
    arg({'repeatable','RepeatableResults'},1,[],'Produce repeatable (vs. randomized) results. If nonzero, the value is taken as the random seed to use for the performed calculation. This way, different numbers (aside from 0) give different repeatable runs.'), ...
    ... % method to evaluate and its arguments
    arg({'args','TrainingArguments'},{},[],'Arguments to training function. Packed in a cell array. Note: if using the default trainer/tester functions, args must at least specify the learning function to be used, optionally followed by arguments to that learning function, e.g. {''lda''} for linear discriminant analysis or {''logreg'',''lambda'',0.1} for logistic regression with ''lambda'',0.1 passed in as user parameters (see ml_train* functions for options).','type','expression'), ...
    arg({'trainer','TrainingFunction'},'@ml_train',[],'Training function. Receives a partition of the data (as produced by the specified partitioner), possibly some further arguments as specified in args, and returns a model (of any kind).','type','expression'), ...
    arg({'tester','PredictionFunction'},'@utl_default_predict',[],'Prediction function. Receives a partition of the data (as produced by the partitioner) and a model (as produced by the trainer), and returns a prediction for every index of the input data, in one of the output formats permitted by ml_predict.','type','expression'), ...
    ... % support for custom data representations
    arg({'target','TargetFunction'},'@utl_default_target',[],'Target extraction function. A function to derive the target variable from a partition of the data (as produced by the partitioner), for evaluation; the allowed format is anything that may be output by ml_predict.','type','expression'), ...
    arg({'partitioner','PartitioningFunction'},'@utl_default_partitioner',[],'Partitioning function. See documentation for the function contract.','type','expression'), ...
    ... % parallel computing settings
    arg({'engine_cv','ParallelEngine','engine'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine to use. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'pool','WorkerPool'},'global',[], 'Worker pool to use. This is typically a cell array, but can also be the string ''gobal'', which stands for the currently globally set up worker pool (see global tracking variable).','type','expression'), ...
    arg({'policy','ReschedulingPolicy'},'global',[], 'Rescheduling policy. This is the name of the rescheduling policy function that controls if and when tasks are being rescheduled. If set to global, the current global setting will be used.'), ...    
    ... % misc arguments
    arg({'cache_fold_results','CacheFoldResults'},false,[],'Whether to cache the per-fold results. This is meant to be used when running very long-running computations on machines that crash frequently enough that partial results need to be saved. In this case, any previously computed results will be loaded from disk.'), ...
    arg({'only_cached_results','OnlyCachedResults'},false,[],'Load only results that are in the cache. This will not run any computations (aside from pre-checks, that can be disabled by setting NoPrechecks to true).'), ...
    arg({'no_prechecks','NoPrechecks'},false,[],'Skip pre-checks that access the data. This can save some time when it would take very long to load the data, especially when performing parallel computation.'), ...
    arg({'collect_models','CollectModels'},false,[],'Collect models per fold. Note that this increases the amount of data returned.'));

data = opts.data; opts = rmfield(opts,'data');

% --- input validation ---

% validate the partitioner argument
if ~isa(opts.partitioner,'function_handle')
    error('The given partitioner must be a function handle, but was: %s',hlp_tostring(opts.partitioner,10000)); end
try
    indexset = opts.partitioner(data,[]);
catch e
    error('The given partitioner failed to calculate the index set (arguments Data,[]) with error: %s',hlp_handleerror(e));
end
if isnumeric(indexset) && isscalar(indexset)
    try
        wholedata = opts.partitioner(data,1:indexset); %#ok<NASGU>
    catch e
        error('The given partitioner failed to partition the data (arguments Data,1:N) with error: %s',hlp_handleerror(e));
    end
elseif iscell(indexset)    
    if isempty(indexset) || any(~cellfun('isclass',indexset,'cell')) || any(cellfun('prodofsize',indexset)<2)
        error('The index-set format returned by the partitioner is unsupported (needs to be either a scalar or a cell array of 2-element cells, but was: %s.',hlp_tostring(indexset,10000)); end
    try
        opts.partitioner(data,indexset{1}{1});
    catch e
        error('The given partitioner failed to partition the data (arguments Data,train-partition-1) with error: %s.',hlp_handleerror(e));
    end
else
    error('The partitioner returned an unsupported index set format: %s',hlp_tostring(indexset,10000));
end

% derive partition index ranges from CV scheme & indexset
inds = make_indices(opts.scheme,indexset,opts.repeatable);

% validate the target argument
if ~isa(opts.target,'function_handle')
    error('The given target argument must be a function handle, but was: %s',hlp_tostring(opts.target,10000)); end
try
    if ~opts.no_prechecks
        tmptargets = opts.target(data); end
catch e
    error('The given target function failed to extract target values from the data with error: %s',hlp_handleerror(e));
end

% parse and validate the metric argument
if ischar(opts.metric) && ~isempty(opts.metric) && opts.metric(1) == '@'
    opts.metric = eval(opts.metric); end
if isempty(opts.metric) || ischar(opts.metric) || (iscell(opts.metric) && all(cellfun(@ischar,opts.metric)))
    opts.metric = make_metric(opts.metric); end
if ~has_stats(opts.metric)
    opts.metric = @(T,P)add_stats(opts.metric(T,P)); end
if ~opts.no_prechecks
    try
        [measure,stats] = opts.metric(tmptargets,tmptargets);
    catch e
        error('Failed to apply the loss metric to target values: %s',hlp_handleerror(e));
    end
    if ~isscalar(measure) && isnumeric(measure)
        error('The given loss metric returned its measure in an unsupported format (should be numeric scalar, but was: %s)',hlp_tostring(measure)); end
    if ~isstruct(stats)
        error('The given loss metric returned its statistics in an unsupported format (should be struct, but was: %s)',hlp_tostring(stats)); end
end
% validate misc arguments
if ~iscell(opts.args)
    error('The given args argument must be a cell array, but was: %s',hlp_tostring(opts.args,10000)); end
if ~isa(opts.trainer,'function_handle')
    error('The given trainer argument must be a function handle, but was: %s',hlp_tostring(opts.trainer,10000)); end
if ~isa(opts.tester,'function_handle')
    error('The given tester argument must be a function handle, but was: %s',hlp_tostring(opts.tester,10000)); end


% --- cross-validation ---

measure = NaN;
stats = struct();    
if ~isempty(inds)
    time0 = tic;
 
    % generate parallelizable tasks for each fold
    for p = length(inds):-1:1
        tasks{p} = {@utl_evaluate_fold,opts,data,inds{p}}; end

    % schedule the tasks
    results = par_schedule(tasks, 'engine',opts.engine_cv, 'pool',opts.pool, 'policy',opts.policy);
    
    % remove empty results (if we're loading in partial results)
    if opts.only_cached_results
        results = results(~cellfun('isempty',results)); end
    
    % collect results
    for p=length(inds):-1:1
        targets{p} = results{p}{1};
        predictions{p} = results{p}{2};
        if opts.collect_models
            models{p} = results{p}{3}; end
    end
        
    % compute aggregate metric / stats
    [dummy,stats] = opts.metric(utl_aggregate_results(targets{:}),utl_aggregate_results(predictions{:})); %#ok<ASGLU>

    % compute per-fold metric / stats
    for p=length(targets):-1:1
        if ~(isempty(targets{p}) && isempty(predictions{p}))
            [measure(p),stats.per_fold(p)] = opts.metric(targets{p},predictions{p}); 
        else
            % this can happen when we are asked to load incomplete results from disk
            fprintf('Note: fold %i had empty predictions and targets.\n',p);
            measure(p) = NaN;
        end
    end
    
    % attach basic summary statistics
    if any(isnan(measure)) && any(~isnan(measure))
        measure = measure(~isnan(measure)); end
    if ~isfield(stats,'measure')
        stats.measure = 'loss'; end
    stats.(stats.measure) = mean(measure);
    stats.([stats.measure '_mu']) = mean(measure);
    stats.([stats.measure '_std']) = std(measure);
    stats.([stats.measure '_med']) = median(measure);
    stats.([stats.measure '_mad']) =  median(abs(measure-median(measure)));
    stats.([stats.measure '_N']) = length(measure);
    measure = mean(measure);
    
    % add per-fold targets & predictions
    for p=1:length(targets)
        stats.per_fold(p).targ = targets{p};
        stats.per_fold(p).pred = predictions{p};
    end

    % add per-fold index sets
    if all(cellfun(@(inds) length(inds{1}) < 10000 && length(inds{2}) < 10000,inds))
        for p=1:length(targets)
            stats.per_fold(p).indices = inds{p}; end
    end                

    % add per-fold models, if any
    if opts.collect_models
        for p=1:length(models)
            stats.per_fold(p).model = models{p}; end
    end
    
    % add additional stats
    stats.time = toc(time0);
end
end


% --- index set generation for data partitioning ---

function inds = make_indices(S,N,repeatable)
% Inds = make_indices(Scheme,Index-Cardinality)
% make cross-validation indices for each fold, from the scheme and the index set cardinality

if isnumeric(N) && isscalar(N)
    % set parameter defaults
    k = 10;                     % foldness or fraction of holdout data
    repeats = 1;                % # of monte carlo repartitions
    randomized = 1;             % randomized indices or not
    margin = 0;                 % width of the index margin between training and evaluation sets
                                % can also be a fraction of the trials if 0<margin<1
    subblocks = 1;              % number of sub-blocks within which to cross-validate

    % parse scheme grammar
    if isnumeric(S)
        % one or more numbers
        switch length(S)
            case 1
                % 0 (skipped CV)
                if S == 0
                    inds = {};
                    return;
                else                    
                    % "folds" format
                    k = S;
                end
            case 2
                % "[repeats, folds]" format
                repeats = S(1);
                k = S(2);
            case 3
                % "[repeats, folds, margin]" format
                repeats = S(1);
                k = S(2);
                margin = S(3);
            otherwise
                error('Unsupported format for cross-validation scheme: %s',hlp_tostring(S));
        end
    elseif ischar(S)
        switch S
            case 'loo'                    
                % "'loo'" format
                k = N;
            case 'trainerr'
                % index set for computing the training-set error
                inds = {{1:N,1:N}};
                return;
            otherwise
                error('Unsupported format for cross-validation scheme: %s',hlp_tostring(S));
        end
    elseif iscell(S) && ~isempty(S)
        if all(cellfun('isclass',S,'cell')) && all(cellfun('prodofsize',S)==2)
            % direct specification of index sets "{{train-inds,test-inds}, {train-inds,test-inds}, ...}"
            inds = S; 
            return;
        elseif ischar(S{1}) && all(cellfun('isclass',S(2:end),'double') | cellfun('isclass',S(2:end),'single')) && all(cellfun('prodofsize',S(2:end))==1)
            % specification as {'string', scalar, scalar, ...}
            switch S{1}
                case {'chron','block'}
                    % "{'chron', k, m}" format
                    randomized = 0;
                    if length(S) > 1
                        k = S{2}; end
                    if length(S) > 2
                        margin = S{3}; end
                    if length(S) == 1 || length(S) > 3
                        error('The {''chron'',...} format must have 1 or 2 numeric parameters; but got: %s',hlp_tostring(S)); end
                case {'subchron','subblock'}
                    % "{'subchron', b, k, m}" format
                    randomized = 0;
                    if length(S) > 1 && isscalar(S{2})
                        subblocks = S{2}; end
                    if length(S) > 2 && isscalar(S{3})
                        k = S{3}; end
                    if length(S) > 3 && isscalar(S{4})
                        margin = S{4}; end
                    if length(S) == 1 || length(S) > 4
                        error('The {''subchron'',...} format must have between 1 and 3 numeric parameters; but got: %s',hlp_tostring(S)); end
                otherwise
                    error('Unsupported format for cross-validation scheme: %s',hlp_tostring(S));
            end
        else
            error('Unsupported format for cross-validation scheme: %s',hlp_tostring(S));
        end
    else
        error('Unsupported cross-validation scheme format: %s',hlp_tostring(S));
    end
    
    if margin > 0 && margin < 1
        margin = round(margin * N); end
    
    % sanity-check parameters
    if k > 1 && round(k) ~= k
        error('The number of folds (k) must be an integer (or a fraction between 0 and 1 to hold out a fraction of the data).'); end
    if k < 0
        error('The number of folds must be nonnegative.'); end
    if round(repeats) ~= repeats
        error('The number of repeats must be an integer.'); end
    if repeats < 1
        error('The number of repeats must be positive.'); end
    if round(margin) ~= margin || margin < 0
        error('The given margin must be a nonnegative integer.'); end
    if round(subblocks) ~= subblocks || subblocks < 1
        error('The given number of sub-blocks for sub-block cross-valdiatio nmust be a positive integer.'); end
        
    % initialize random number generation
    if randomized && repeatable
        if hlp_matlab_version < 707
            % save & override RNG state
            randstate = rand('state'); %#ok<*RAND>
            rand('state',5182+repeatable);
        else
            % create a legacy-compatible RandStream
            randstream = RandStream('swb2712','Seed',5182+repeatable);
        end
    end

    % generate evaluation index sets from the parameters
    try
        % set up permutation
        if subblocks == 1
            perm = 1:N;
        elseif subblocks <= N
            Nblock = N + subblocks-mod(N,subblocks);
            perm = reshape(1:Nblock,[],subblocks)';
            perm = round(perm(:)*N/Nblock);
        else
            error('There are more subblocks than trials in the data.');
        end

        % build indices
        inds = {};
        for r=1:repeats
            if randomized
                if hlp_matlab_version < 707
                    perm = randperm(N);
                else
                    perm = randstream.randperm(N);
                end
            end
            if k < 1
                % p-holdout
                inds{end+1} = sort(perm((round(N*(1-k))+1):N)); %#ok<AGROW>
            else
                % k-fold
                for i=0:k-1
                    inds{end+1} = sort(perm(1+floor(i*N/k) : min(N,floor((i+1)*N/k)))); end %#ok<AGROW>
            end
        end
    catch err
        % % this error is thrown only after the subsequent delicate RNG state restoration
        indexgen_error = err;
    end
    
    % restore saved RNG state
    if randomized && repeatable && hlp_matlab_version < 707        
        rand('state',randstate); end;
    
    % optionally throw any error that happened during previous index set creation
    if exist('indexgen_error','var')
        rethrow(indexgen_error); end 

    % add complementary training index sets, handling the margin
    for i=1:length(inds)
        tmpinds = true(1,N);
        tmpinds(inds{i}) = 0;
        for j=1:margin
            tmpinds(max(1,inds{i}-j)) = 0;
            tmpinds(min(N,inds{i}+j)) = 0;
        end
        inds{i} = {find(tmpinds),inds{i}}; %#ok<AGROW>
    end
elseif iscell(N)
    % the partitioner returned the index sets already; pass them through
    inds = N;
else
    error('Unsupported index set format: %s',hlp_tostring(N));
end
end


% test whether the given metric supplies stats or not
function result = has_stats(metric)
try 
    [x,y] = metric([],[]);  %#ok<NASGU,ASGLU>
    result = true;
catch e
    result = ~any(strcmp(e.identifier,{'MATLAB:TooManyOutputs','MATLAB:maxlhs','MATLAB:unassignedOutputs'}));
end
end



% add stats to the result of some metric
function [result,stats] = add_stats(result)
stats.value = result;
end

% create a cross-validation metric
% (this is a sub-function to keep the created anonymous function lean to serialize)
function func = make_metric(metric)
    func = @(T,P) ml_calcloss(metric,T,P);
end
