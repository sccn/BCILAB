function id = onl_newpredictor(name, model, streams, predict_at)
% Create a new predictor from a predictive model, and tie it to some stream(s).
% Id = onl_newpredictor(PredictorName,Model,StreamNames,PredictAt)
%
% A predictor is created from a predictive model (which was previously calibrated via bci_train or
% the GUI), and linked to streams, from which it reads data. Once data is available (after some data
% was supplied via onl_append()), the predictor can be queried for its prediction, given the most
% recently supplied data.
% 
% In:
%   PredictorName : name of the predictor to be created; a variable of this name will be created in 
%                   the MATLAB base workspace to hold the predictor's state.
%
%   Model : A model data structure (as obtained from bci_train) based on which the predictor shall be 
%           created; typically this is a model struct, but for convenience it can be a file name, 
%           variable name in the base workspace, or a cell array of {file name, variable name} to 
%           refer to a variable inside a .mat file. The model is not modified by this function.
%
%   SourceStreamNames : Optional names of stream data structures in the MATLAB base workspace to
%                       consider as possible data sources (previously created with onl_newstream); 
%                       if a stream contains all channels that are needed by the predictor, or 
%                       alternatively has the right number and type of channels it will be considered 
%                       as a potential source stream unless ambiguous.
%
%   PredictAt : Determines where predictions are made when onl_predict is called; if {}, the
%               prediction is made at the most recently added sample of the stream, and if nonempty
%               a prediction is made for every marker in the stream that matches the given cell
%               array of event types (which may also contain wildcard characters).
%
% Out:
%   Id : a unique id number for the predictor; same as PredictorName.predictorid
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
    % for each pipeline...
    for c=1:length(predictor.tracking.filter_graph)
        % substitute the node set_makepos('signal',x,...,'online_epoching',y,...)
        % by set_makepos('Signal',set_targetmarkers(x,predict_at),...,'online_epoching','at_targets',...)
        [match,pos] = utl_find_filter(predictor.tracking.filter_graph{c},'set_makepos');
        % override the online_epoching argument
        argpos = find(strcmp(match.parts,'online_epoching'))+1;
        if isempty(argpos)
            error('Your model was apparently calculated with a BCILAB distribution that predated the marker streaming capability; please recompute it.'); end
        match.parts{argpos} = 'at_targets';
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
