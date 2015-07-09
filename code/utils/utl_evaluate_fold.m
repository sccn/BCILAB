function result = utl_evaluate_fold(opts,data,inds)
% Internal to utl_crossval; Learns on a training set and tests on a test set.
% The training set is indexed by inds{1}, and the test set is indexed by inds{2}.
%
% In:
%   Options : options struct as used for utl_crossval
%
%   Data    : a representation of some (indexable) data
%
%   Indices : a cell array of {training indices, test indices}
%             note: optionally, a third element may be supplied which is taken as a cell array of 
%             additional arguments to the trainer
%
% Out:
%   Result  : a cell array of {real target values, predicted target values}
%
% See also:
%   utl_crossval
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-07

% utl_evaluate_fold_version<1.0> -- for the cache
dp;

disp('fold.');
if opts.cache_fold_results
    result = hlp_diskcache({'cvfolds' 'load_only', opts.only_cached_results}, ...
        @cached_evaluate,rmfield(opts,'only_cached_results'),data,inds);
    % return empty result record if we're in cache loadonly mode
    if strcmp(result,'hlp_diskcache:notfound')   
        result = []; end
else
    result = cached_evaluate(opts,data,inds);
end


function result = cached_evaluate(opts,data,inds)
% this function is just an exception-catching wrapper to evaluate_internal
% cached_evaluate_version<1.0> -- for the cache
if ~opts.tolerate_exceptions
    result = evaluate_internal(opts,data,inds);
else
    try
        result = evaluate_internal(opts,data,inds);    
    catch e
        fprintf('utl_evaluate_fold: suppressing exception: %s\n',hlp_handleerror(e));
        result = [];
    end
end
    

function result = evaluate_internal(opts,data,inds)
% the function that does the actual work

trainset = opts.partitioner(data,inds{1});

% learn a model on the training partition
if length(inds) >= 3
    model = opts.trainer(trainset,opts.args{:},inds{3}{:});
else
    model = opts.trainer(trainset,opts.args{:});
end

% and record both the known targets (result{1}) and predictions (result{2}) on the test partition,
% and optionally the trained model (result{3})
testset = opts.partitioner(data,inds{2});
result = {opts.target(testset), opts.tester(testset,model), quickif(opts.collect_models,model,struct())};
