function run_readdataset(varargin)
% Receive (simulated) real-time data from a dataset.
% run_readdataset(MatlabStream,Dataset,UpdateFrequency)
%
% In:
%   StreamName : Name of the stream; a variable with this name will be created in the MATLAB workspace 
%                to hold the stream's data. If such a variable already exists it will be overridden.
%
%   Dataset : Dataset to play back. This is the EEGLAB dataset that shall be played back in real
%             time. Can also be the name of a variable in the MATLAB workspace. (default: 'lastdata')
%
%   UpdateFrequency : The rate at which new chunks of data is polled from the device, in Hz. 
%                     (default: 20)
%
%   BufferLength : Internal buffering length. This is the maximum amount of backlog that you 
%                  can get. (default: 10)
%
%   ConvertToDouble : Convert to double. Always convert the signal to double precision.
%                     (default: true)
%
% Examples:
%   % open a new input stream, and update it with data read in real time from an EEGLAB data set
%   run_readdataset('mystream',EEG)
%
%   % as before, but pass the arguments by name
%   run_readdataset('MatlabStream','mystream','Dataset',EEG)
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','Dataset');

% read arguments...
arg_define(varargin, ...
    arg({'new_stream','StreamName','MatlabStream'}, 'laststream',[],'MATLAB Stream Name. A variable with this name will be created in the MATLAB workspace to hold the stream''s data. If such a variable already exists it will be overridden.'), ...
    arg({'in_dataset','Dataset'},'lastdata',[],'Dataset to play back. This is the EEGLAB dataset that shall be played back in real time. Can also be the name of a variable in the MATLAB workspace.','type','expression'), ...
    arg({'update_freq','UpdateFrequency'},20,[0 Inf],'Update frequency. The rate at which new chunks of data is polled from the device, in Hz.'), ...
    arg({'buffer_len','BufferLength'},10,[],'Internal buffering length. This is the maximum amount of backlog that you can get.'), ...
    arg({'always_double','ConvertToDouble'},true,[],'Convert to double. Always convert the signal to double precision.'), ...
    arg({'looped_playback','LoopedPlayback','Looping','looping'},true,[],'Looped playback. Whether the playback is looping (otherwise the stream will stop at some point).'));

% evaluate dataset, if it's an expression
in_dataset = exp_eval_optimized(in_dataset);
in_dataset.starttime = tic;
if ~isempty(in_dataset.event)
    in_dataset.event_latencies = round([in_dataset.event.latency]);
else
    in_dataset.event_latencies = [];
end
    

% open a new online stream
onl_newstream(new_stream,rmfield(in_dataset,'data'),'buffer_len',buffer_len);

% start a background reading job
onl_read_background(new_stream,@(stream)read_block(in_dataset,stream,always_double,looped_playback), update_freq);

disp('Now reading...');


% background block reader function
function result = read_block(in_dataset,stream,always_double,looped_playback)
% get current position and read-out range
curpos = round(toc(in_dataset.starttime)*in_dataset.srate);
if curpos > size(in_dataset.data,2) && ~looped_playback
    error('BCILAB:EndOfStream','End of stream.'); end    
range = 1+mod(stream.smax : curpos-1,size(in_dataset.data,2));
% get the next data block
block = in_dataset.data(:,range);
if always_double
    block = double(block); end
% get events
if ~isempty(range) && ~isempty(in_dataset.event_latencies)
    events = in_dataset.event(in_dataset.event_latencies>=range(1)&in_dataset.event_latencies<=range(end));
    if ~isempty(events)
        [events.latency] = arraydeal([events.latency]-range(1)+1);end
else
    events = [];
end
result = {block,events};
