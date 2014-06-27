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


% learn a model
trainset = opts.partitioner(data,inds{1});
if length(inds) >= 3
    model = opts.trainer(trainset,opts.args{:},inds{3}{:});
else
    model = opts.trainer(trainset,opts.args{:});
end

% and record both the real test targets (result{1}) and the prediction on the test data (result{2})
testset = opts.partitioner(data,inds{2});
result = {opts.target(testset), opts.tester(testset,model),quickif(opts.collect_models,model,struct())};