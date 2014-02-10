function [measure,stats] = utl_crossval(data, varargin)
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
%   Data : some data that can be partitioned using index sets, such as, for example,
%          {X,y} with X being the [NxF] training data and y being the [Nx1] labels 
%          (N=#trials, F=#features)
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
%                          formats that can be produced by ml_predict (default: @ml_predict)
%
%               'scheme': cross-validation scheme, can be one of the following formats (default: 10)
%                * 0: skip CV, return NaN and empty statistics
%                * k: k-fold randomized CV (or, if 0<k<1, k-holdhout CV)
%                * [r k]: r times repartitioned k-fold randomized CV (or, if 0<k<1, k-holdout CV 
%                         (with k a fraction))
%                * [r k m]: r times repartitioned k-fold randomized CV (or, if 0<k<1, k-holdout CV 
%                           (with k a fraction)) with m indices margin width (between training and
%                           test set)
%                * 'loo': leave-one-out CV
%                * {'chron', k} or {'block', k}: k-fold chronological/blockwise CV
%                * {'chron', k, m} or {'block', k, m}: k-fold chronological/blockwise CV with m 
%                                                      indices margin width (between training and 
%                                                      test set)
%                * 'trainerr': This is the training-set error, i.e. it is a measure of the 
%                              separability of the data (not of the generalization ability of the 
%                              classifier); it is usually an error to report this number in a paper
%
%               'partitioner': partitioning function for the data, receives two parameters: 
%                              (data, index vector)
%                              * if the index vector is empty, should return the highest index in 
%                                the data OR a cell array of {training-partition,test-partition}
%                                for each fold (in the latter case, the partitioner fully controls
%                                the cross-validation scheme); for each fold, these two outputs
%                                will be fed into the partitioner as second argument to generate the 
%                                training set and the test set, respectively
%                              * otherwise, it should return data subindexed by the index vector
%                              default: provides support for cell arrays, numeric arrays, struct
%                                       arrays and {Data,Target} cell arrays 
%
%               'target': a function to derive the target variable from a partition of the data 
%                         (as produced by the partitioner), for evaluation; the allowed format is
%                         anything that may be output by ml_predict default: provides support for
%                         {Data,Target} cell arrays
%
%               'metric': metric to be employed, applied both to results of each fold and results 
%                         aggregated over all folds
%                         * function handle: a custom, user-supplied loss function; receives target 
%                           data in the first argument and prediction data in the second argument;
%                           each can be in any format that can be produced by ml_predict (but can be
%                           expected to be mutually consistent). shall return a real number
%                           indicating the summary metric over all data, and optionally additional
%                           statistics in a struct
%                         * string: use ml_calcloss, with 'metric' determining the loss type
%                         * default/empty: use 'mcr','mse','nll','kld' depending on supplied target 
%                           and prediction data formats, via ml_calcloss
%
%               'args': optional arguments to the training function, packed in a cell array 
%                       (default: empty)
%                       note: if using the default trainer/tester combination, args must at least 
%                             specify the learning function to be used, e.g. {'lda'} for linear
%                             discriminant analysis (see ml_train* functions for options)
%
%               'repeatable': whether the randomization procedure shall give repeatable results 
%                             (default: 1); different numbers (aside from 0) give different
%                             repeatable runs, i.e. the value determines the randseed
%
%               'return_models': whether to return models trained for each fold (default: false)
%
%               'engine_cv': parallelization engine to use (default: 'global'); see par_beginsschedule
%                   
%               'pool'  : worker pool to use (default: 'global'); see par_beginsschedule
%
%               'policy': scheduling policy to use (default: 'global'); see par_beginschedule
%
% Out:
%   Measure : a measure of the overall performance of the trainer/tester combination, w.r.t. to the 
%             target variable returned by the target function computed by the metric
%
%   Stats   : additional statistics, as produced by the metric
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

% read arguments
opts = hlp_varargin2struct(varargin, ...
    'scheme',10, ...
    'partitioner',@utl_default_partitioner,  ...
    'target',@utl_default_target, ...
    'metric',@(T,P) ml_calcloss([],T,P),  ...
    'trainer',@ml_train,  ...
    'tester',@utl_default_predict,  ...
    'args',{}, ...
    'repeatable',1,  ...
    'collect_models',false, ...
    'engine_cv','global', ...
    'pool','global', ...
    'policy','global');

% string arguments are considered to be variants of the default metric
if isempty(opts.metric) || ischar(opts.metric) || (iscell(opts.metric) && all(cellfun(@ischar,opts.metric)))
    opts.metric = @(T,P)ml_calcloss(opts.metric,T,P); end
if ~has_stats(opts.metric)
    opts.metric = @(T,P)add_stats(opts.metric(T,P)); end

% derive indices from CV scheme & N
inds = make_indices(opts.scheme,opts.partitioner(data,[]),opts.repeatable);

if isempty(inds)
    measure = NaN;
    stats = struct();
else
    
    time0 = tic;
 
    % generate tasks for each fold
    for p = 1:length(inds)
        tasks{p} = {@utl_evaluate_fold,opts,data,inds{p}}; end

    % schedule the tasks
    results = par_schedule(tasks, 'engine',opts.engine_cv, 'pool',opts.pool, 'policy',opts.policy);
    
    % collect results
    for p=1:length(inds)
        % collect results
        targets{p} = results{p}{1};
        predictions{p} = results{p}{2};
        if opts.collect_models
            models{p} = results{p}{3}; end
    end
        
    % compute aggregate metric
    [dummy,stats] = opts.metric(utl_aggregate_results(targets{:}),utl_aggregate_results(predictions{:})); %#ok<ASGLU>

    % add per-fold metric    
    for p=1:length(targets)
        stats.per_fold(p) = hlp_getresult(2,opts.metric,targets{p},predictions{p}); end
    
    % calculate basic cross-validation statistics
    tmp = [stats.per_fold.(stats.measure)];
    stats.(stats.measure) = mean(tmp);
    stats.([stats.measure '_mu']) = mean(tmp);
    stats.([stats.measure '_std']) = std(tmp);
    stats.([stats.measure '_med']) = median(tmp);
    stats.([stats.measure '_mad']) =  median(abs(tmp-median(tmp)));
    measure = mean(tmp);
    
    % also add the original targets & predictions
    for p=1:length(targets)
        stats.per_fold(p).targ = targets{p};
        stats.per_fold(p).pred = predictions{p};
    end

    % add original index sets
    if all(cellfun(@(inds) length(inds{1}) < 10000 && length(inds{2}) < 10000,inds))
        for p=1:length(targets)
            stats.per_fold(p).indices = inds{p}; end
    end                

    % add collected models, if any
    if opts.collect_models
        for p=1:length(models)
            stats.per_fold(p).model = models{p}; end
    end
    
    % add additional stats
    stats.time = toc(time0);
    
end
end



function inds = make_indices(S,N,repeatable)
% Inds = make_indices(Scheme,Index-Cardinality)
% make cross-validation indices for each fold, from the scheme and the index set cardinality

if iscell(N) && all(cellfun('isclass',N,'cell')) && all(cellfun('prodofsize',N)>=2)
    % the partitioner returned the index sets already; pass them through
    inds = N;
else
    if strcmp(S,'trainerr')
        % special case: index set for computing the training-set error
        inds = {{1:N,1:N}};
    else
        % regular case: proper cross-validation
        
        % set parameter defaults
        k = 10;                     % foldness or fraction of holdout data
        repeats = 1;                % # of monte carlo repartitions
        randomized = 1;             % randomized indices or not
        margin = 0;                 % width of the index margin between training and evaluation sets
        subblocks = 1;              % number of sub-blocks within which to cross-validate
        
        % parse scheme grammar
        if isnumeric(S) && length(S) == 3
            % "[repeats, folds, margin]" format
            repeats = S(1);
            k = S(2);
            margin = S(3);
        elseif isnumeric(S) && length(S) == 2
            % "[repeats, folds]" format
            repeats = S(1);
            k = S(2);
        elseif iscell(S) && ~isempty(S) && ischar(S{1}) && ...
                (strcmp('chron',S{1}) || strcmp('block',S{1}))
            % "{'chron', k, m}" format
            randomized = 0;
            if length(S) > 1 && isscalar(S{2})
                k = S{2}; end
            if length(S) > 2 && isscalar(S{3})
                margin = S{3}; end
        elseif iscell(S) && ~isempty(S) && ischar(S{1}) && ...
                (strcmp('subchron',S{1}) || strcmp('subblock',S{1}))
            % "{'subchron', b, k, m}" format
            randomized = 0;
            if length(S) > 1 && isscalar(S{2})
                subblocks = S{2}; end
            if length(S) > 2 && isscalar(S{3})
                k = S{3}; end
            if length(S) > 3 && isscalar(S{4})
                margin = S{4}; end
        elseif isscalar(S)
            % "k" format
            k = S;
        elseif ischar(S) && strcmp(S,'loo')
            % "'loo'" format
            k = N;
        elseif iscell(S) && all(cellfun('isclass',S,'cell')) && all(cellfun('prodofsize',S)>=2)
            % direct specification
            inds = S; 
            return;
        else
            error('unknown cross-validation scheme format');
        end
        
        % check for skipped CV
        inds = {};
        if k <= 0
            return; end
        
        if randomized && repeatable
            if hlp_matlab_version < 707
                % save & override RNG state
                randstate = rand('state'); %#ok<*RAND>
                rand('state',5182+repeatable); %#ok<RAND>
            else
                % create a legacy-compatible RandStream
                randstream = RandStream('swb2712','Seed',5182+repeatable);
            end
        end
        
        % generate evaluation index sets from the parameters
        try
            if subblocks == 1
                perm = 1:N;
            else
                Nblock = N + subblocks-mod(N,subblocks);
                perm = reshape(1:Nblock,[],subblocks)';
                perm = round(perm(:)*N/Nblock);
            end
            
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
                    inds{end+1} = {sort(perm(1+(0:(round(N*k)-1))))};
                else
                    % k-fold
                    for i=0:k-1
                        inds{end+1} = sort(perm(1+floor(i*N/k) : min(N,floor((i+1)*N/k)))); end
                end
            end
        catch err
            % % this error is thrown only after the subsequent delicate RNG state restoration
            indexgen_error = err;
        end
        
        if randomized && repeatable && hlp_matlab_version < 707
            % restore saved RNG state
            rand('state',randstate); %#ok<RAND>
        end
        
        if exist('indexgen_error','var')
            rethrow(indexgen_error); end % throw the error that happened during previous index set creation
        
        % add complementary training index sets
        for i=1:length(inds)
            tmpinds = true(1,N);
            tmpinds(inds{i}) = 0;
            for j=1:margin
                tmpinds(max(1,inds{i}-j)) = 0;
                tmpinds(min(N,inds{i}+j)) = 0;
            end
            inds{i} = {find(tmpinds),inds{i}};
        end
    end
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

