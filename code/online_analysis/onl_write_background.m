function onl_write_background(varargin)
% Periodically process data using a predictive model, and write results to some external device.
% onl_write_background(ResultWriter,StreamNames,Model,OutputFormat,UpdateFrequency,StartDelay,PredictorName,PredictAt,Verbose,EmptyResultValue)
% 
% This is a convenience function which simplifies the definition of components which load and then
% periodically query a predictive model, in real time, and forward the results to some external
% device. The function is internally implemented using a timer that periodically triggers the
% computation of updated estimates, and their transfer to the data sink.
%
% In:
%   ResultWriter : Callback function that receives one or more BCI estimates and writes them to some
%                  external device. The format passed to the ResultWriter is according to OutputFormat.
%
%   StreamNames : optional names of stream data structures in the MATLAB base workspace to consider as
%                 possible data sources (previously created with onl_newstream); any stream that 
%                 contains channels that are needed by the predictor will be linked to it (assuming 
%                 that the choice of stream to use is not ambiguous). 
%
%                 The identification of needed channels is primarily done on the basis of the channel
%                 labels -- if a stream has channels with labels that are required by a filter pipeline,
%                 it will be used as a source for this pipeline. The framework attempts to gracefully
%                 handle cases where a stream only provides a subset of the channels that were in the 
%                 training set and the model only effectively operates on this subset via flt_selchans.
%
%   Model : A model data structure (as obtained from bci_train) according to which predictions shall be 
%           made; typically this is a model struct, but for convenience it can be a file name, 
%           variable name in the base workspace, or a cell array of {file name, variable name} to 
%           refer to a variable inside a .mat file. The model is not modified by this function. 
%           (default: 'lastmodel')
%
%   OutputFormat : Output data format, see onl_predict (default: 'distribution')
%
%   UpdateFrequency : Frequency at which the device should be queried, in Hz (default: 10)
%
%   StartDelay : Delay before real-time processing begins; grace period until user resources are 
%                created (default: 1)
%
%   PredictorName : Name of a predictor data structure that shall be created in the MATLAB base 
%                   workspace to hold the predictor's state. If a variable of this name already exists
%                   it will be overridden. (default: 'lastpredictor')
%
%   PredictAt : Predict at markers. If nonempty, this is a cell array of online target markers 
%               relative to which predictions shall be made. If empty, predictions are always made 
%               on the most recently added sample. (default: {})
%
%   Verbose : Verbose output. If false, the console output of the online pipeline will be suppressed.
%             (default: false)
%
%   EmptyResultValue : Empty-result value. This value is returned for predictions that yielded no 
%                      result (e.g., due to an error or because not enough data was available).
%                      (default: NaN)
%
% Examples:
%   % after a predictive model has been learned using bci_train, and a data stream supplying raw
%   % data has been established, load the model into the online system and periodically send its 
%   % outputs to a target destination
%   onl_write_background(@send_outputs_to_destination,'mystream')
%
%   % as before, but also specify a custom output format and a higher update frequency
%   onl_write_background(@send_outputs_to_destination,'mystream','lastmodel','expectation',25)
%
%   % as before, but pass all arguments by their short names
%   onl_write_background('ResultWriter',@send_outputs_to_destination,'MatlabStream','mystream','Model','lastmodel','OutputFormat','expectation','UpdateFrequency',25)
%
% See also:
%   onl_predict
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-18

% read options
arg_define(varargin, ...
    arg_norep({'result_writer','ResultWriter'},[],[],'Result writing callback. Callback function that receives one or more BCI estimates and writes them to some external device. The format passed to the ResultWriter is according to OutputFormat.'), ...
    arg({'in_stream','StreamName','MatlabStream'}, 'laststream',[],'Input stream name. This is the name of the stream data structure in the MATLAB base workspace that shall be read from. Can also be a cell array of stream names, if multiple, or empty if non-ambiguous.','typecheck',false,'shapecheck',false), ...
    arg({'pred_model','Model'}, 'lastmodel', [], 'Predictive model. A model data structure (as obtained from bci_train) according to which predictions shall be made; typically this is a model struct, but for convenience it can be a file name, variable name in the base workspace, or a cell array of {file name, variable name} to refer to a variable inside a .mat file. The model is not modified by this function.','type','expression'), ...
    arg({'out_form','OutputFormat'},'distribution',{'expectation','distribution','mode','raw'},'Format of the produced output values. Can be the expected value (posterior mean) of the target variable, or the distribution over possible target values (probabilities for each outcome, or parametric distribution), or the mode (most likely value) of the target variable.'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. This is the rate at which the outputs should be calculated.'), ...
    arg({'start_delay','StartDelay'}, 1, [0 Inf],'Start-up delay. Delay before real-time processing begins; grace period to initialize everything.'), ...
    arg({'pred_name','PredictorName'}, 'lastpredictor',[],'Name of new predictor. Name of a predictor data structure that shall be created in the MATLAB base workspace to hold the predictor''s state. If a variable of this name already exists it will be overridden.'), ...
    arg({'predict_at','PredictAt'}, {},[],'Predict at markers. If nonempty, this is a cell array of online target markers relative to which predictions shall be made. If empty, predictions are always made on the most recently added sample.','type','expression'), ...
    arg({'verbose','Verbose'}, false,[],'Verbose output. If false, the console output of the online pipeline will be suppressed.'), ...
    arg({'empty_result_value','EmptyResultValue'},NaN,[],'Empty-result value. This value is returned for predictions that yielded no result (e.g., due to an error or because not enough data was available).','type','expression'));

% input validation
if ischar(result_writer) && ~isempty(result_writer) %#ok<NODEF>
    if result_writer(1) ~= '@' && ~exist(result_writer,'file')
        error('The given ResultWriter argument (%s) is not a valid function name',result_writer); end
    result_writer = str2func(result_writer); 
end
if ~isa(result_writer,'function_handle')
    error('The given ResultWriter argument must be a function handle (or function name), but was: %s',hlp_tostring(result_writer,10000)); end
if ~iscell(in_stream) %#ok<NODEF>
    in_stream = {in_stream}; end
for c=1:length(in_stream)
    stream_name = in_stream{c};
    if ~isvarname(stream_name)
        error('The given StreamName argument must be a valid variable name, but was: %s',hlp_tostring(stream_name,10000)); end
    try
        stream = evalin('base',stream_name);
    catch e
        error('Failed to look up stream named %s in MATLAB base workspace with error: %s',stream_name,e.message);
    end
    if ~isstruct(stream)
        error('The given data structure named %s in the MATLAB base workspace was expected to be a stream data structure, but was not a struct (wrong name?): %s',stream_name,hlp_tostring(stream,10000)); end
    if ~isfield(stream,'streamid')
        if isfield(stream,{'data','srate'})
            error('The given stream data structure named %s appears to be an EEGLAB data set struct but is not a stream (use onl_newstream to create a valid stream)',stream_name);
        else
            error('The given data structure named %s is not a valid stream (use onl_newstream to create a valid stream)',stream_name);
        end
    end
    streamids{c} = stream.streamid; %#ok<AGROW>
end

% create new predictor
predid = onl_newpredictor(pred_name,pred_model,in_stream,predict_at);

% create & start timer (which periodically writes to the stream)
start(timer('ExecutionMode','fixedRate', 'Name',[pred_name '_timer'], 'Period',1/update_freq, ...
    'StartDelay',start_delay, 'TimerFcn',@(timer_handle,varargin) write_data(pred_name,in_stream,out_form,result_writer,predid,streamids,timer_handle,verbose,empty_result_value)));

% background data writer
function write_data(predictor,stream_names,fmt,result_writer,pred_id,stream_ids,timer_handle,verbose,empty_result_value)
try
    % check if the stream and the predictor are still there
    for c=1:length(stream_names)
        s = evalin('base',stream_names{c});
        if s.streamid ~= stream_ids{c}
            error('Note: the stream named %s was recreated.',stream_names{c}); end
    end
    p = evalin('base',predictor);
    if p.predictorid ~= pred_id
        error('Note: the predictor named %s was recreated.',predictor); end
    % make a prediction
    y = onl_predict(predictor,fmt,~verbose,empty_result_value);
    % and write it out
    try
        result_writer(y);
    catch e
        disp('Error in result-writing function:');
        hlp_handleerror(e);
    end
catch e
    if ~strcmp(e.identifier,'MATLAB:UndefinedFunction')
        hlp_handleerror(e); end    
    % stream or predictor have changed (e.g., replaced/deleted) --> stop timer
    stop(timer_handle);
    delete(timer_handle);
end
