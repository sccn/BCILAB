function stream_id = onl_read_background(varargin)
% Read from an external device in the background and periodically update a stream with the results.
% onl_read_background(MatlabStream, BlockReader, UpdateFrequency)
%
% This is a convenience function which simplifies the implementation of a data source that runs in 
% parallel to the data processing done by predictive models. It is internally based on a timer which 
% periodically invokes a user-supplied data reading function.
% 
% In:
%   StreamName : Name of the stream data structure in the MATLAB workspace that shall be updated 
%                in the background (must have been previously created via onl_newstream) 
%                (default: 'laststream')
% 
%   BlockReader : Callback function that reads a block from the device (format: [Channels x Samples])
%                 if no data is available, this function should return an empty result. It may also
%                 return a cell array of the form {[Channels x Samples],TimeStamp} or {[Channels x
%                 Samples],Markers,TimeStamp}. More generally, all returned cells are used as
%                 arguments to onl_predict. Optionally, this function may take the current stream
%                 variable as input.
%
%   UpdateFrequency : Frequency at which the device should be queried, in Hz (default: 10)
%
% Example:
%   % after a stream has been openend, ...
%   onl_newstream('mystream','srate',200,'chanlocs',{'C3','Cz','C4'});
%
%   % ensure that it gets updated periodically (here: at 30 Hz) using new samples from some device
%   % the device function should return all samples that have accumulated since it was last called.
%   onl_read_background('mystream', @get_samples_from_my_device,30);
%
% See also:
%   onl_append, onl_newstream
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-01-18

% read options
arg_define(varargin, ...
    arg({'stream_name','StreamName','MatlabStream'}, 'laststream',[],'Stream name to create. Name of the stream data structure in the MATLAB workspace that shall be updated in the background (must have been previously created via onl_newstream).'), ...
    arg({'block_reader','BlockReader'},[],[],'Block-reading function. Callback function that reads a block from the device (format: [Channels x Samples]) if no data is available, this function should return an empty result. It may also return a cell array of the form {[Channels x Samples],TimeStamp} or {[Channels x Samples],Markers,TimeStamp}. More generally, all returned cells are used as arguments to onl_predict. Optionally, this function may take the current stream variable as input.','type','expression'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. New data is polled at this rate, in Hz.'));

% input validation
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
stream_id = stream.streamid;
if ischar(block_reader) && ~isempty(block_reader)
    if block_reader(1) ~= '@' && ~exist(block_reader,'file')
        error('The given BlockReader argument (%s) is not a valid function name',block_reader); end
    block_reader = str2func(block_reader); 
end
if ~isa(block_reader,'function_handle')
    error('The given BlockReader argument must be a function handle (or function name), but was: %s',hlp_tostring(block_reader,10000)); end


% create & start timer
start(timer('ExecutionMode','fixedRate', 'Name',[stream_name '_timer'], 'Period',1/update_freq, ...
    'TimerFcn',@(timer_handle,varargin)append_data(stream_name,stream_id,timer_handle,block_reader)));

% timer callback: append data to a stream
function append_data(stream_name,stream_id,timer_handle,read_block)
try
    % check if the stream is still there
    x = evalin('base',stream_name);
    if x.streamid ~= stream_id
        error('Note: the stream named %s was recreated.',stream_name); end
    try
        % get a new block
        if nargin(read_block) == 1
            block = read_block(x);
        else
            block = read_block();
        end
        % append it to the stream
        if iscell(block)        
            onl_append(stream_name,block{:});
        else
            onl_append(stream_name,block);
        end
    catch e
        if strcmp(e.identifier,'BCILAB:EndOfStream')
            disp('Encountered end-of-stream.');
            evalin('base',['clear ' stream_name]);
            stop(timer_handle);
            delete(timer_handle);
        else
            disp('Error in block-reading function:');
            hlp_handleerror(e);
        end
    end
catch e
    if ~strcmp(e.identifier,'MATLAB:UndefinedFunction')
        hlp_handleerror(e); end
    % stream has been changed (e.g., replaced/deleted) --> stop timer
    stop(timer_handle);
    delete(timer_handle);
end
