function model = utl_complete_model(model,func,bundle)
% Internal. Complete the definition of a predicitive model with some meta-data.
% utl_complete_model(Model, PredictionFunction, Bundle)
%
% This function ensures that model structs created by BCI paradigm plugins have a uniform format
% and contain the set of fields needed by downstream processing functions (e.g., the offline and 
% online prediction framework). This function relieves the BCI paradigms from specifying each and
% every of those fields manually, by automatically appending fields that can be deduced from the
% data or other parts of the model.
%
% This function also performs some early validation of the model structure to catch possible errors
% early on rather than deep into a subsequent processing run.
%
% In:
%   Model : a predictive model that will be finalized by this function; may have the following fields:
%           .tracking.filter_graph : cell array of symbolic filter expressions, one per inlet of the
%                                    prediction function; can be used to filter streams or bundles
%
%                                    Can also be passed either as an EEGLAB data set (with .tracking
%                                    field) or as a stream bundle; in this case, the processing that
%                                    has been applied to those data sets will be used as a template
%                                    to lay out the filter graph.
%
%                                    If not passed, will be filled in by the Bundle argument.
%
%           .tracking.prediction_function : function that takes a window of the output of filter graph's
%                                           output and computes a prediction from it; If not passed,
%                                           will be filled in based on the PredictionFunction argument 
%
%           .tracking.prediction_channels : cell array of chanlocs structures expected by the prediction
%                                           function; has one element for each inlet of the prediction 
%                                           function.
%
%                                           If not passed, will be filled in based on the Bundle argument.
%
%
%   PredictionFunction : handle of the prediction function (supplied by the framework code that calls 
%                        utl_complete_model)
%
%   Bundle : a stream bundle fully preprocessed by the filter graph, used to extract some
%            meta-data (supplied by framework code)
%
% Out:
%   Model : a completed model struct; the following fields are ensured to exist:
%           .tracking                --> tracking information about the model for offline/online processing
%            .filter_graph           --> a cell array of (symbolic) filter expressions, denoting the 
%                                        filter operations that yield each stream inlet that is 
%                                        accepted by the model's prediction function. Each expression
%                                        may contain leaf nodes of the type @rawdata, which are free
%                                        inlets of the filter graph (tagged with the channels expected
%                                        at that inlet)
%            .prediction_function    --> a function that takes a stream bundle with one stream for 
%                                        each entry of the filter graph; produces a prediction for each
%                                        epoch/segment in the data (therefore works both online and offline)
%            .prediction_channels    --> a cell array of chanlocs structs that defines that channels
%                                        expected by the prediction function at each of its inlets
%                                        (this is redundant additional information)
%
% See also:
%   bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-28
dp;

% input validation
if ~isstruct(model) || ~isscalar(model)
    error('The given model must be a 1x1 struct, but was: %s',hlp_tostring(model,10000)); end
if exist('bundle','var')
    if ~isstruct(bundle) || ~isscalar(bundle)
        error('The given stream bundle must be a 1x1 struct, but was: %s',hlp_tostring(bundle,10000)); end
    if ~isfield(bundle,'streams') && isfield(bundle,'tracking') && isfield(bundle.tracking,'online_expression')
        bundle.streams = {bundle}; end
    for s=1:length(bundle.streams)
        utl_check_fields(bundle.streams{s},{'tracking','chanlocs'},'signal','signal'); end
end
if exist('func','var') && ~(isa(func,'function_handle') || isa(func,'char'))
    error('The given PredictionFunction argument must either be a function handle or a string.'); end

% the .tracking field holds the data that is of interest to the offline/online evaluation system
% it is allowed for the paradigm to fill in the following field manually, in order to yield custom behavior
if ~isfield(model,'tracking')
    model.tracking = struct(); end

% add the filter_graph, if necessary
if ~isfield(model.tracking,'filter_graph')
    if ~exist('bundle','var')
        error('Cannot automatically deduce the .tracking.filter_graph field for this model, because the data on which its state shall be initialized is not known or available;\nthis most likely the case if your model was calibrated on a dataset collection and the framework could not decide which dataset in the collection to initialize the filters on.\nYour can set the field directly in your paradigm''s calibrate() function.'); end
    for s=1:length(bundle.streams)
        model.tracking.filter_graph{s} = bundle.streams{s}.tracking.online_expression; end
else
    % a filter graph was passed in; check format
    if all(isfield(model.tracking.filter_graph,{'tracking','chanlocs'}))
        % ... as a data set; override the bundle with it
        bundle = struct('streams',{{model.tracking.filter_graph}});        
        % set the filter graph accordingly
        model.tracking.filter_graph = {bundle.streams{1}.tracking.online_expression};
    elseif isfield(model.tracking.filter_graph,'streams')
        % ... as a stream bundle: override the bundle with it
        bundle = model.tracking.filter_graph;
        % and initialize the filter graph accordingly
        model.tracking.filter_graph = {};
        for s=1:length(bundle.streams)
            model.tracking.filter_graph{s} = bundle.streams{s}.tracking.online_expression; end
    elseif all(isfield(model.tracking.filter_graph,{'head','parts'}))
        % ... as an abstract expression
        model.tracking.filter_graph = {model.tracking.filter_graph};
    elseif iscell(model.tracking.filter_graph)
        % ... in the proper format (cell array)
        if isempty(model.tracking.filter_graph)
            error('The given model''s .tracking.filter_graph is empty.'); end
        for c=1:length(model.tracking.filter_graph)
            if ~all(isfield(model.tracking.filter_graph{c},{'head','parts'}))
                error('The given model''s .tracking.filter_graph is not a valid filter graph (at least one of its cells is lacking the required fields .head and/or .parts): %s',hlp_tostring(model.tracking.filter_graph)); end
        end
    else
        error('The given model''s .tracking.filter_graph is not a valid filter graph: %s',hlp_tostring(model.tracking.filter_graph));
    end
end

% add the prediction function, if necessary
if ~isfield(model.tracking,'prediction_function')
    if ~exist('func','var')
        error('If the given model is lacking a .tracking.prediction_function field, the prediction function should be supplied as second input to utl_complete_model.'); end
    model.tracking.prediction_function = func; 
end

% validate the model's prediction function
if ischar(model.tracking.prediction_function)
    if strncmp(model.tracking.prediction_function,'Paradigm',8)
        try
            % Paradigm class reference: try to instantiate
            instance = eval(model.tracking.prediction_function); %#ok<NASGU>
        catch e
            error('Failed to instantiate the paradigm referred to by the model''s .tracking.prediction_function with error: %s',hlp_handleerror(e));
        end
    else
        if ~exist(model.tracking.prediction_function,'file')
            error('The function referred to by the model''s .tracking.prediction_function field does not exist: %s',model.tracking.prediction_function); end
    end
else
    tmp = char(model.tracking.prediction_function);
    if tmp(1) ~= '@' && ~exist(tmp,'file')
        error('The function referred to by the model''s .tracking.prediction_function field does not exist: %s',char(model.tracking.prediction_function)); end
end

% add and validate the channel locations (these are technically redundant, as they can be 
% derived from the preprocessing chain) (deduced from the bundle)
try
    if ~isfield(model.tracking,'prediction_channels')
        if ~exist('bundle','var')
            error('Cannot automatically assign the .tracking.prediction_channels field because the data set from which these channels shall be determined could not be deduced by the framework; consider setting this field manually in the calibrate() function of your paradigm.'); end
        for s=1:length(bundle.streams)
            if ~((isstruct(bundle.streams{s}.chanlocs) && isfield(bundle.streams{s}.chanlocs,'labels')) || iscellstr(bundle.streams{s}.chanlocs))
                error('The .chanlocs fields in one of the streams of the stream bundle is malformed; needs to be a chanlocs struct (or cell array of channel labels), but was: %s',s,hlp_tostring(bundle.streams{s}.chanlocs)); end
            if iscellstr(bundle.streams{s}.chanlocs)
                bundle.streams{s}.chanlocs = struct('labels',bundle.streams{s}.chanlocs); end
            model.tracking.prediction_channels{s} = bundle.streams{s}.chanlocs; 
        end
    else
        for s=1:length(model.tracking.prediction_channels)
            if ~((isstruct(model.tracking.prediction_channels{s}) && isfield(model.tracking.prediction_channels{s},'labels')) || iscellstr(model.tracking.prediction_channels{s}))
                error('The given model''s .tracking.prediction_channels{k} field is malformed for k=%i; needs to be a chanlocs struct or cell array of channel labels, but was: %s',s,hlp_tostring(model.tracking.prediction_channels{s})); end 
        end
    end
catch e %#ok<CTCH>
    fprintf('Could not initialize prediction channels: %s\n',hlp_handleerror(e));
end

if length(model.tracking.filter_graph) ~= length(model.tracking.prediction_channels)
    error('The given filter graph has a different number of elements than the given prediction channels; they should have the same length.'); end

% finally add also a timestamp (so that we can sort them by creation time)
model.timestamp = now;
