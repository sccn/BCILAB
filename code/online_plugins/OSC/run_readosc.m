function run_readosc(varargin)
% Receive real-time data from OSC.
% run_readosc(MatlabStream,InputPort,InputPath,InputMetadata,UpdateFrequency)
%
% Note: The address that is used to pass data to this interface must be identical to InputPath
%       (i.e. there is no pattern matching implemented).
%
% In:
%   StreamName : Name of the stream; a variable with this name will be created in the MATLAB workspace 
%                to hold the stream's data. If such a variable already exists it will be overridden.
%
%   InputPort : Port at which to listen for OSC input data (default: 12345)
%
%   InputPath : OSC path at which to listen for input data (default: '/data')
%
%   InputMetadata : Meta-data of the input stream. This is a struct or cell array of name-value
%                   pairs with meta-data fields to use. The mandatory fields are 'srate' and
%                   'chanlocs', where chanlocs is either a channel locations struct array, or
%                   a cell array of channel names, or the number of channels (in which case
%                   a cell array of the form {'A1','A2', ..., 'A32','B1','B2', ...} is created).
%                   Optionally, the field 'datasource' can be set to point to a dataset on disk or
%                   to a MATLAB workspace variable, to serve as the source of meta-data.
%
%   UpdateFrequency : The rate at which new chunks of data is polled from the device, in Hz. 
%                     (default: 20)
%
% Example:
%   % open an OSC input stream, listening at port 22050 under address /data
%   run_readosc('mystream',22050,{'srate',512,'chanlocs','mychanlocs.sfp'})
%
%   % as before, but pass the arguments by name
%   run_readosc('MatlabStream','mystream','InputPort',22050,'InputMetadata',{'srate',512,'chanlocs','mychanlocs.sfp'})
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','OSC');

% read options
arg_define(varargin, ...
    arg({'new_stream','StreamName','MatlabStream'}, 'laststream',[],'MATLAB Stream Name. A variable with this name will be created in the MATLAB workspace to hold the stream''s data. If such a variable already exists it will be overridden.'), ...
    arg({'in_port','InputPort'}, 12345,uint32([1 65535]),'Input OSC listen port. This is the network port at which the stream will listen for data.'), ...
    arg({'in_path','InputPath'}, '/data',[],'Input OSC method path. This is the method/address name at which to listen for input data.'), ...
    arg_sub({'in_metadata','InputMetadata'},{},@utl_parse_metadata, 'Meta-data of the input stream. These are fields as they appear in EEGLAB data sets; only sampling rate and channel labels are mandatory.'), ...
    arg({'update_freq','UpdateFrequency'},25,[0 Inf],'Update frequency. The rate at which new chunks of data is polled from the device, in Hz.'));

% parse the meta-data spec
meta = utl_parse_metadata(in_metadata);

% open a DataRiver input stream (and make sure it gets cleaned up eventually)
if ~exist('osc_new_server','file')
    try
        build_osc;
    catch
        error('The OSC library has not been built for your platform yet; see dependencies/OSC* for more info.'); 
    end
end
server = osc_new_server(in_port);
deleter = onCleanup(@()osc_free_server(server));

% create online stream stream (with all metadata from the dataset)
onl_newstream(new_stream,meta);

% create background reading job
onl_read_background(new_stream,@()read_block(server,deleter,in_path),update_freq)

% background data reading function
function block = read_block(server,deleter,ipath)
messages = osc_recv(server,0);
buffer = {};
for m=1:length(messages)
    if strcmp(messages{m}.path,ipath)
        buffer{m} = [messages{m}.data{:}]'; end
end
block = [buffer{:}];
