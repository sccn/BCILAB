function stream_id = onl_read_background(varargin)
% Read from an external device in the background and periodically update a stream with the results.
% onl_read_background(MatlabStream, BlockReader, UpdateFrequency)
%
% This is a convenience function which simplifies the implementation of a data source that runs in 
% parallel to the data processing done by predictive models. It is internally based on a timer which 
% periodically invokes a user-supplied data reading function.
% 
% In:
%   MatlabStream : name of the stream to be updated (default: 'laststream')
%                  (must have been previously created via onl_newstream)
% 
%   BlockReader : function that reads a block from the device (format: [Channels x Samples])
%                 if no data is available, this function should return an empty result.
%                 It may also return a cell array of the form {[Channels x Samples],TimeStamp}. More 
%                 generally, all returned cells are used as arguments to onl_predict.
%                 Optionally, this function may take the current stream variable as input.
%
%   UpdateFrequency : frequency at which the device should be queried, in Hz (default: 25)
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
    arg({'stream_name','MatlabStream'}, 'laststream',[],'New Stream to create. This is the name of the stream within the MATLAB environment.'), ...
    arg({'block_reader','BlockReader'},'randn(32,1)',[],'Block-reading function. This function reads a new block from some device','type','expression'), ...
    arg({'update_freq','UpdateFrequency'},10,[],'Update frequency. New data is polled at this rate, in Hz.'));

% get the stream's id
stream_id = evalin('base',[stream_name '.streamid']);

% create & start timer
start(timer('ExecutionMode','fixedRate', 'Name',[stream_name '_timer'], 'Period',1/update_freq, ...
    'TimerFcn',@(timer_handle,varargin)append_data(stream_name,stream_id,timer_handle,block_reader)));

% timer callback: append data to a stream
function append_data(stream_name,stream_id,timer_handle,read_block)
try
    % check if the stream is still there
    x = evalin('base',stream_name);
    if x.streamid ~= stream_id
        error('Stream changed.'); end
    % get a new block
    if nargin(read_block) == 1
        block = read_block(x);
    else
        block = read_block();
    end
    % append it to the stream
    if iscell(block)        
        % includes a timestamp
        onl_append(stream_name,block{:});
    else
        % raw data block
        onl_append(stream_name,block);
    end
catch e
    if ~strcmp(e.identifier,'MATLAB:UndefinedFunction')
        env_handleerror(e); end
    % stream has been changed (e.g., replaced/deleted) --> stop timer
    stop(timer_handle);
    delete(timer_handle);
end
