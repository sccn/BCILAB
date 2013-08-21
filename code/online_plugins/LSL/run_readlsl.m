function run_readlsl(varargin)
% Receive real-time data from a source on the lab streaming layer.
% run_readlsl(MatlabStream,SelectionProperty,SelectionValue,UpdateFrequency)
%
% This plugin connects to and receives data from a device on the lab streaming layer. The device
% is identified by one of its properties (e.g., name or type).
%
% In:
%   MatlabStream : name of the stream to create in the MATLAB environment (default: 'laststream')
%
%   SelectionProperty : property that shall be used to find the desired device (default: 'type')
%
%   SelectionValue : the value of the chosen property that the device must have (default: 'EEG')
%
%   ConvertToDouble : Always convert the signal to double precision. (default: true)
%   
%   UpdateFrequency : this is rate at which new data is polled from the device, in Hz (default: 100)
%
%   BufferLength : Internal buffering length. This is the maximum amount of backlog that you can
%                  get, in seconds. (default: 30)
%
%   ChannelOverride : Override channel labels. This allows to replace the channel labels that 
%                     are provided by the stream. (default: {})
%
% Examples:
%   % receive data from a device that contains gaze data
%   run_readlsl('mystream','type','Gaze');
%
%   % read from an EEG stream (default) but use custom channel labels
%   run_readlsl('ChannelOverride', {'C3','C4','Cz','O1','O2'};
%
%   % read from a stream that has the name 'PhaseSpace'
%   run_readlsl('SelectionProperty','name','SelectionValue','PhaseSpace')
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2012-03-21

persistent lib;

% declare the name of this component (shown in the menu)
declare_properties('name','Lab streaming layer');

% read options
opts = arg_define(varargin, ...
    arg({'new_stream','MatlabStream'}, 'laststream',[],'New Stream to create. This is the name of the stream within the MATLAB environment.'), ...
    arg({'property','SelectionProperty'}, 'type',[],'Selection property. The selection criterion by which the desired device is identified on the net. This is a property that the desired device must have (e.g., name, type, desc/manufacturer, etc.'), ...
    arg({'value','SelectionValue'}, 'EEG',[],'Selection value. This is the value that the desired device must have for the selected property (e.g., EEG if searching by type, or Biosemi if searching by manufacturer).'), ...
    arg({'always_double','ConvertToDouble'},true,[],'Convert to double. Always convert the signal to double precision.'), ...
    arg({'update_freq','UpdateFrequency'},20,[],'Update frequency. New data is polled at this rate, in Hz.'), ...
    arg({'buffer_len','BufferLength'},30,[],'Internal buffering length. This is the maximum amount of backlog that you can get.'), ...
    arg({'channel_override','ChannelOverride'}, [], [], 'Override channel labels. This allows to replace the channel labels that are provided by the stream.','type','cellstr','shape','row'));

% get a library handle (here with an explicit path because we want it to work if the toolbox is compiled, too)
if isempty(lib)
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin')); end

% look for the desired device
disp(['Looking for a device with ' opts.property '=' opts.value ' ...']);
result = {};
while isempty(result)
    result = lsl_resolve_byprop(lib,opts.property,opts.value); end

% create a new inlet
disp('Opening an inlet...');
inlet = lsl_inlet(result{1});

% query the stream info...
info = inlet.info();

% check srate
if info.nominal_srate() <= 0
    warning('BCILAB:run_readlsl:not_supported','The given stream has a variable sampling rate. This plugin currently only supports data with regular sampling rate.'); end

if isempty(opts.channel_override)
    % try to get the channel labels & check them
    channels = {};    
    ch = info.desc().child('channels').child('channel');
    if ch.empty()
        ch = info.desc().child('channel'); end
    while ~ch.empty()
        name = ch.child_value_n('label');
        if isempty(name)
            name = ch.child_value_n('name'); end
        if name
            channels{end+1} = name; end
        ch = ch.next_sibling_n('channel');
    end
else
    channels = opts.channel_override;
end

if length(channels) ~= info.channel_count()
    disp('The number of channels in the steam does not match the number of labeled channel records. Using numbered labels.');
    channels = cellfun(@(k)['Ch' num2str(k)],num2cell(1:info.channel_count(),1),'UniformOutput',false);
end

    
% create online stream (using appropriate meta-data)
srate = info.nominal_srate();
if srate == 0
    % for non-uniformly sampled streams we assume a sampling rate that is approx in the ballpark of 
    % 100 Hz in order to determine the buffer length
    srate = 100; 
end 
onl_newstream(opts.new_stream,'srate',srate,'chanlocs',channels,'buffer_len',opts.buffer_len);

% start background acquisition
if opts.always_double
    onl_read_background(opts.new_stream,@()double(inlet.pull_chunk()),opts.update_freq);
else
    onl_read_background(opts.new_stream,@()inlet.pull_chunk(),opts.update_freq);
end
