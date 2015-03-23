function outbundle = utl_preprocess_bundle(inbundle,model)
% Internal. Pre-process a stream bundle with a model's filter graph.
%
% BCI models have a "filter graph" (possibly trivial), which represents a directed acyclical graph
% of filter stages with their parameters. This function takes a multi-stream signal (a stream
% bundle) and applies the filter graph to it, yielding a new multi-stream bundle that's the output
% of the graph.
%
% In:
%   Bundle : a stream bundle to pre-process (struct with a field .streams that's a cell array of 
%            EEGLAB dataset structs)
%
%   Model : predictive model to use
%
% Out:
%   Bundle : processed stream bundle
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-11-23
dp;

% input validation
if ~isstruct(model) || ~isscalar(model)
    error('The given model argument must be a 1x1 struct.'); end
if ~isfield(model,'tracking')
    error('The given model data structure is lacking the required .tracking field; its fields are: %s',hlp_tostring(fieldnames(model))); end
if ~isfield(model.tracking,'filter_graph')
    error('The given model data structure is lacking the required .tracking.filter_graph field.'); end
if ~iscell(model.tracking.filter_graph) || ~all(cellfun(@(f)all(isfield(f,{'head','parts'})),model.tracking.filter_graph))
    error('The given model''s .tracking.filter_graph has an unsupported structure (should be a cell array of expressions, but was: %s',hlp_tostring(model.tracking.filter_graph,10000)); end
if ~isfield(model.tracking,'prediction_channels')
    error('The given model data structure is lacking the required .tracking.prediction_channels field.'); end
if ~iscell(model.tracking.prediction_channels) || ~all(cellfun(@(c)isfield(c,'labels'),model.tracking.prediction_channels))
    error('The given model''s .tracking.prediction_channels has an unsupported structure (should be a cell array of chanlocs structs, but was: %s',hlp_tostring(model.tracking.filter_graph,10000)); end
if length(model.tracking.filter_graph) ~= length(model.tracking.prediction_channels)
    error('The given model''s .tracking.filter_graph and .tracking.prediction_channels have different lengths. This is not permitted.'); end

% first resolve all rawdata nodes in the model's filter graph according to the input bundle
resolved_graph = utl_resolve_streams(model.tracking.filter_graph,inbundle,model.tracking.prediction_channels);

% then evaluate each chain in the graph and store the results as an output stream bundle
outbundle.streams = cell(1,length(resolved_graph));
for g=1:length(resolved_graph)
    outbundle.streams{g} = exp_eval_optimized(resolved_graph{g}); end
