function predictions = utl_collection_tester(testcollection,model,predict_func)
% Internal. Apply a predictive model to a dataset collection, as part of a cross-validation.
%
% See also:
%   bci_train, utl_collection_partition, utl_collection_targets

for k=1:length(testcollection)
    predictions{k} = predict_func(utl_preprocess_bundle(testcollection{k},model),model); end
predictions = utl_aggregate_results(predictions{:});
