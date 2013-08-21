function outbundle = utl_preprocess_bundle(inbundle,model)
% Internal. Pre-process a stream bundle with a model's filter graph.

% first resolve all @rawdata nodes in the model's filter graph according to the input bundle
resolved_graph = utl_resolve_streams(model.tracking.filter_graph,inbundle,model.tracking.prediction_channels);
% then evaluate each chain in the graph and store the results as an output stream bundle
outbundle.streams = cell(1,length(resolved_graph));
for g=1:length(resolved_graph)
    outbundle.streams{g} = exp_eval_optimized(resolved_graph{g}); end
