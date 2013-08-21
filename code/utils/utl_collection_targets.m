function targets = utl_collection_targets(testcollection)
% Internal. Obtain the target values for a dataset collection, as part of a cross-validation.
%
% See also:
%   utl_collection_tester

% note: these are drawn from the first stream (as all streams have the same markers)
for k=1:length(testcollection)
    targets{k} = set_gettarget(testcollection{k}.streams{1}); end
targets = utl_aggregate_results(targets{:});
