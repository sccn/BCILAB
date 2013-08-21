function model = utl_complete_model(model,func,bundle)
% Internal. Complete the definition of a predicitive model with some meta-data.
% utl_complete_model(Model, PredictionFunction, Bundle)
% In:
%   Model : a predictive model that will be finalized by this function; may have the following fields:
%           .tracking.filter_graph : cell array of symbolic filter expressions, one per inlet of the
%                                    prediction function; can be used to filter streams or bundles
%
%                                    Can also be passed either as an EEGLAB data set (with .tracking
%                                    field) or as a stream bundle; in this case, the processing that
%                                    has been applied to those data sets will be used as a template
%                                    to lay out to set up the filter graph.
%
%                                    If not passed, will be filled in by the Bundle argument.
%
%           .tracking.prediction_function : function that takes a window of the output of filter graph's
%                                           output and computes a prediction from it; If not passed,
%                                           will be filled in based on the PredictionFunction argument 
%
%           .tracking.prediction_window : window length of data expected by the prediction_function, 
%                                         in samples (if 0, the most recent chunk of data will be passed);
%                                         has one element for each inlet of the prediction function;
%
%                                         If not passed, will be filled in based on the Bundle argument.
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
%            .prediction_window      --> an array that determines, for each inlet of the prediction function, 
%                                        the number of samples expected by it in the epochs/segments,
%                                        or 0 if the chunk of most recent data in the stream is expected.
%            .prediction_channels    --> a cell array of chanlocs structs that defines that channels
%                                        expected by the prediction function at each of its inlets
%                                        (this is redundant additional information)
%
% See also:
%   bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-28

% the .tracking field holds the data that is of interest to the offline/online evaluation system
% it is allowed for the paradigm to fill in the following field manually, in order to yield custom behavior
if ~isfield(model,'tracking')
    model.tracking = struct(); end

% add the filter_graph, if necessary
if ~isfield(model.tracking,'filter_graph')
    if ~exist('bundle','var')
        error('Cannot automatically deduce the .tracking.filter_graph field for this model, because the data on which its state shall be initialized is not known or available; this most likely the case if your model was calibrated on a dataset collection and the framework could not decide which dataset in the collection to initialize the filters on. Your may set the field directly in your paradigm''s calibrate() function.'); end
    for s=1:length(bundle.streams)
        model.tracking.filter_graph{s} = bundle.streams{s}.tracking.online_expression; end
else
    % a filter graph was passed in; check format
    if all(isfield(model.tracking.filter_graph,{'data','chanlocs'}))
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
    else
        % in the regular format
        if ~iscell(model.tracking.filter_graph)
            model.tracking.filter_graph = {model.tracking.filter_graph}; end
    end
end

% add the prediction function, if necessary
if ~isfield(model.tracking,'prediction_function')
    model.tracking.prediction_function = func; end

try
    % add the prediction_window, if necessary (deduced from the bundle)
    if ~isfield(model.tracking,'prediction_window')
        if ~isempty(bundle.streams{1}.epoch) || size(bundle.streams{1}.data,3) > 1
            % data is epoched: by default the window length is the length of the epochs of the respective streams in the bundle
            for s=1:length(bundle.streams)
                model.tracking.prediction_window(s) = size(bundle.streams{1}.data,2); end
        else
            % data is not epoched: by default the window length is the length of the most recent chunk in each stream
            model.tracking.prediction_window = zeros(1,length(bundle.streams));
        end
    end
    
    % add the channel locations (these are technically redundant, as they can be derived from the preprocessing chain)
    % (deduced from the bundle)
    if ~isfield(model.tracking,'prediction_channels')
        if ~exist('bundle','var')
            errore('Cannot automatically assign the .tracking.prediction_channels field because the data set from which these channels shall be determined could not be deduced by the framework; consider setting this field manually in the calibrate() function of your paradigm. This is a one-time warning.'); end
        for s=1:length(bundle.streams)
            model.tracking.prediction_channels{s} = bundle.streams{s}.chanlocs; end
    end
catch
end

% finally add also a timestamp (so that we can sort them by creation time)
model.timestamp = now;
