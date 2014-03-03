function run_readbiosemi(varargin)
% Receive real-time data from BioSemi.
% run_readbiosemi(MatlabStream,ChannelRange,SamplingRate,UpdateFrequency)
%
% This plugin connects to a BioSemi ActiveTwo amplifier (Mk1 or Mk2), using the BioSemi USB driver.
% It was tested on Linux and Window, but it may be necessary to recompile the driver interface for
% your platform (esp. on Linux); see the readme files in dependencies/BioSemi-2010-11-19 for this.
%
% The meta-data (channel names and order) is pre-defined by the amplifier, but it is possible to 
% read only a subset of the provided data (for efficiency), using the ChannelRange and SamplingRate
% parameters.
%
% In:
%   StreamName : Name of the stream; a variable with this name will be created in the MATLAB workspace 
%                to hold the stream's data. If such a variable already exists it will be overridden.
%
%   ChannelRange : Numeric vector of channel indices that should be recorded (referring to the 
%                  default BioSemi channel order); default: 3:131
%
%   SamplingRate : Sampling rate for the amplifier, in Hz (default: 256)
%
%   BufferLength : Internal buffering length. This is the maximum amount of backlog that you can
%                  get, in seconds. (default: 30)
%
%   UpdateFrequency : The rate at which new chunks of data is polled from the device, in Hz. 
%                     (default: 10)
%
% Examples:
%   % open a biosemi input stream that is sampled at 512 Hz and updated at 30 Hz 
%   % (using the default name for the stream)
%   run_readbiosemi('UpdateFrequency',30, 'SamplingRate',512);
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','BioSemi amplifier');

% read options
opts = arg_define(varargin, ...
    arg({'new_stream','StreamName','MatlabStream'}, 'laststream',[],'MATLAB Stream Name. A variable with this name will be created in the MATLAB workspace to hold the stream''s data. If such a variable already exists it will be overridden.'), ...
    arg({'channel_range','ChannelRange'}, 3:128+3,uint32([1 1000]),'Reduced channel range. Allows to specify a sub-range of the default BioSemi channels.'), ...
    arg({'sample_rate','SamplingRate'}, 256,[1 100000],'Sampling rate. In Hz.'), ...
    arg({'update_freq','UpdateFrequency'},10,[0 Inf],'Update frequency. The rate at which new chunks of data is polled from the device, in Hz.'), ...
    arg({'buffer_len','BufferLength'},10,[],'Internal buffering length. This is the maximum amount of backlog that you can get.'));

% open a BioSemi connection (using an explicit path because it should also work if the toolbox is compiled)
conn = bs_open(env_translatepath('dependencies:/BioSemi-2010-11-19'));

% create online stream (using appropriate meta-data from the connection)
onl_newstream(opts.new_stream,'srate',opts.sample_rate,'chanlocs',conn.channels(opts.channel_range),'xmin',toc(uint64(0)),'buffer_len',opts.buffer_len);

% start background acquisition
onl_read_background(opts.new_stream,@()read_block(conn,opts),opts.update_freq);

disp('Now reading...');



% background block reader function
function block = read_block(conn,opts)
% get a new block from biosemi
block = bs_read(conn);
% decimate it (note: this is a poor implementation!)
if ~isempty(block)
    block = block(opts.channel_range,1:double(conn.srate)/opts.sample_rate:end); end
