function id = onl_newpredictor(name, model, streams)
% Create a new predictor from a predictive model, and tie it to some stream(s).
% Id = onl_newpredictor(Name, Model, Streams)
%
% A predictor is created from a predictive model (which was previously calibrated via bci_train or
% the GUI), and linked to streams, from which it reads data. Once data is available (after some data
% was supplied via onl_append()), the predictor can be queried for its prediction, given the most
% recently supplied data.
% 
% In:
%   Name : name of the predictor to be created; a variable of this name will be created in the base
%          workspace to hold the predictor's data; calling 'clear' in the workspace will erase the
%          predictor.
%
%   Model : A model to be turned into a predictor; this can be a model struct, or a base workspace 
%           variable name, or a file name, or a cell array of {file name, variable name} to refer
%           to a variable inside a file. Models are calibrated via bci_train or the GUI.
%
%   Streams : optional names of streams (previously created with onl_newstream) to consider as
%             possible data sources; any stream that contains channels that are needed by the
%             predictor will be linked to it (assuming that the choice of stream to use is not
%             ambiguous). 
%
%             The identification of needed channels is primarily done on the basis of the channel
%             labels -- if a stream has channels with labels that are required by a filter pipeline,
%             it will be used as a source for this pipeline. The framework attempts to gracefully
%             handle cases where a stream only provides a subset of the channels that were in the 
%             training set and the model only effectively operates on this subset via flt_selchans.
%
% Out:
%   Id : a unique id number for the predictor; same as name.predictorid
%
% Notes:
%   When a detector gets linked to a stream, it will stop to function after the stream has been 
%   deleted. 
%   
%   Appends the following fields to the model: 
%    .pipelines --> a cell array of filter chains (with buffers), each formatted as a recursive 
%                   filter expression 
%    .name --> name of the predictor
%
%   Resolves the rawdata(...) expression into the form rawdata(stream_name) where 
%   stream_name is the name of the stream that is the source.
% 
% See also:
%   onl_newstream, onl_append, onl_predict
%
% Examples:
%   % assuming that the required channels are present in some previously created online stream,
%   % and that a model has been computed previously, create a new predictor from it for online processing
%   onl_newpredictor('mypredictor',lastmodel)
%
%   % if there are multiple applicable streams, directly pass the set of data streams that should be 
%   % used as source data
%   onl_newpredictor('mypredictor',lastmodel,{'mystream1','mystream2'})
%
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-06-21

if ~exist('name','var')
    error('Please specify a name for the predictor.'); end
if ~isvarname(name)
    error('The name of the predictor must be a valid variable name.'); end
if ~exist('model','var')
    error('Please specify a model for the predictor.'); end

% parse the model
if ischar(model)
    % try to evaluate in the workspace
    try
        model = evalin('base',model);
    catch
        % if that fails, try to load it as a file name
        try
            model = io_load(model);
            if ~isfield(model,'tracking') || ~isfield(model.tracking,'prediction_function')
                % the loaded model is lacking the appropriate fields; check if there are variables
                % in the loaded data which are valid models
                candidates = {};
                for f = fieldnames(model)'
                    fname = f{1};
                    if isfield(model.(fname),'tracking') && isfield(model.(fname).tracking,'prediction_function')
                        candidates{end+1} = fname; end %#ok<AGROW>
                end
                if length(candidates) > 1
                    error('BCILAB:onl_newpredictor:ambiguous',['The file given as the model contains multiple candiate variables:\n' ...
                        hlp_tostring(candidates) '; please pass a file or model structure which is non-ambiguous.']); 
                elseif isempty(candidates)
                    error('BCILAB:onl_newpredictor:load_error','The given file contains no valid model.');
                else
                    model = model.(candidates{1});
                end
            end
        catch
            error('BCILAB:onl_newpredictor:load_error','The given model string could not be interpreted (neither as a file name nor as a workspace variable).');
        end
    end
elseif iscell(model) && length(model) == 2 && iscellstr(model)
    % two-element cell-string arrays are interpreted as {filename varname}.
    try
        model = getfield(io_load(model{1},model{2}),model{2}); %#ok<GFLD>
    catch
        error('BCILAB:onl_loaddetector:load_error',['The file ' model{1} ' and/or its variable ' model{2} ' could not be loaded.']);
    end
elseif ~isstruct(model) || isempty(model)
    error('BCILAB:onl_newpredictor:invalid_model','The given data is not a valid model.');
elseif ~isfield(model,'tracking') || ~isfield(model.tracking,'prediction_function') || ~isfield(model.tracking,'filter_graph') || ~isfield(model.tracking,'prediction_window')
    error('BCILAB:onl_newpredictor:invalid_model','The given data structure is not a valid model (lacking required fields).');
end


% get the streams
if ~exist('streams','var') || isempty(streams)
    % find all admissible streams in the workspace....
    vars = evalin('base','whos');
    vars = vars(strcmp({vars.class},'struct'));
    streams = {vars(cellfun(@(x)all(isfield(evalin('base',x),{'buffer','smax'})),{vars.name})).name};
end

if ~iscell(streams)
    streams = {streams}; end

for s=1:length(streams)
    if ~ischar(streams{s})
        error('BCILAB:onl_newpredictor:invalid_streams','The Streams argument must be passed as the names under which the streams were loaded, instead of as structs.'); end
    if ~isvarname(streams{s})
        error('BCILAB:onl_newpredictor:invalid_streams','One of the supplied stream names is not a valid matlab variable name (and thus cannot refer to a stream): %s.',streams{s}); end
end

try
    % create the predictor
    predictor = model;
    predictor.name = name;
    predictor.pipelines = predictor.tracking.filter_graph;
    % resolve the rawdata nodes into the correct stream
    predictor.pipelines = utl_resolve_streams(predictor.pipelines,streams,predictor.tracking.prediction_channels);
    % assign a unique id
    id = fresh_id('bcilab_predictors');
    predictor.predictorid = id;
    % assign to base workspace
    assignin('base',name,predictor);
catch e
    env_handleerror(e);
    error('BCILAB:onl_newpredictor:unexpected','The given model has an unexpected structure.');
end
