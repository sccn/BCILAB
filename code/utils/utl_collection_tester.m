function predictions = utl_collection_tester(testcollection,model,predict_func)
% Internal. Apply a predictive model to a dataset collection, as part of a cross-validation.
% Predictions = utl_collection_tester(TestCollection,Model,PredictionFunction)
%
% In:
%   TestCollection : dataset collection (cell array) to which a model shall be applied; the elements
%                    can be stream bundles or EEGLAB dataset structs
%
%   Model : predictive model to use
%
%   PredictionFunction : prediction function to use
%
% Out:
%   Predictions : the predictions of the model on all recordings, concatenated
%
% See also:
%   bci_train, utl_collection_partition, utl_collection_targets
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-29
dp;

% input validation
if ~iscell(testcollection) || ~all(cellfun('isclass',testcollection,'struct'))
    error('The given TestCollection argument must be a cell array of structs, but was: %s',hlp_tostring(testcollection,10000)); end
if ~isfield(model,'tracking') && all(isfield(model.tracking,{'filter_graph','prediction_channels'}))
    error('The given Model argument must be a predictive model (have field .tracking.filter_graph and .tracking.prediction_channels), but was: %s',hlp_tostring(model,10000)); end
if ~isa(predict_func,'function_handle')
    error('The give PredictionFunction argument must be a function handle.'); end

for k=1:length(testcollection)
    predictions{k} = predict_func(utl_preprocess_bundle(testcollection{k},model),model); end %#ok<AGROW>
predictions = utl_aggregate_results(predictions{:});
