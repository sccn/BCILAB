function id = onl_newpredictor(name, model, streams, predict_at)
% Create a new predictor from a predictive model, and tie it to some stream(s).
% Id = onl_newpredictor(Name,Model,Streams,PredictAt)
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
%   PredictAt : Determines where predictions are made when onl_predict is called; if {}, the
%               prediction is made at the most recently added sample of the stream, and if nonempty
%               a prediction is made for every marker in the stream that matches the given cell
%               array of event types (which may also contain wildcard characters).
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
if ~exist('streams','var')
    streams = []; end
if ~exist('predict_at','var') || isempty(predict_at)
    predict_at = {}; end

% initialize the predictor
predictor = utl_loadmodel(model);
predictor.name = name;

% optionally set up the filter graph to generate epochs relative to desired target markers
if ~isempty(predict_at)
    if ~iscellstr(predict_at)
        error('PredictAt must be a cell array of strings.'); end
    % for each pipeline...
    for c=1:length(predictor.tracking.filter_graph)
        % substitute the node set_makepos('signal',x,...,'online_epoching',y,...)
        % by set_makepos('Signal',set_targetmarkers(x,predict_at),...,'online_epoching','at_targets',...)
        [match,pos] = utl_find_filter(predictor.tracking.filter_graph{c},'set_makepos');
        % override the online_epoching argument
        match.parts{find(strcmp(match.parts,'online_epoching'))+1} = 'at_targets';
        % insert set_targetmarkers stage
        idx = find(strcmp(match.parts,'signal'))+1;
        match.parts{idx} = set_targetmarkers('signal',match.parts{idx},'eventmap',predict_at,'epoch_bounds',[0 0],'eventfield','type','prune_nontarget',false,'avoid_boundaries',false,'arg_direct',true);
        % write back
        predictor.tracking.filter_graph{c} = subsasgn(predictor.tracking.filter_graph{c},pos,match);
    end
end

% instantiate pipelines for each output of the filter graph
for k=1:length(predictor.tracking.filter_graph)
    predictor.pipelines{k} = onl_newpipeline(predictor.tracking.filter_graph{k},streams,predictor.tracking.prediction_channels{k}); end

% convert string-valued prediction functions to a callable function
if ischar(predictor.tracking.prediction_function)
    if strncmp(predictor.tracking.prediction_function,'Paradigm',8)
        % Paradigm class reference: instantiate and get the paradigm's prediction function
        instance = eval(predictor.tracking.prediction_function); %#ok<NASGU>
        predictor.tracking.prediction_function = eval('@instance.predict');
    else
        % cast to a function
        predictor.tracking.prediction_function = str2func(predictor.tracking.prediction_function);
    end
end

% determine whether the prediction function is stateful
predictor.stateful = is_stateful(predictor.tracking.prediction_function,[],[]);

% add an .arg_direct field set to true in order to ensure that, if the prediction
% function uses arg_define, it does not attempt to use slow parsing on its inputs
predictor.arg_direct = 1;

% assign a unique id
id = fresh_id('bcilab_predictors');
predictor.predictorid = id;

% assign to base workspace
assignin('base',name,predictor);
