function run_readdatariver(varargin)
% Receive real-time data from DataRiver.
% run_readdatariver(MatlabStream,DiskStream,InputMetadata,UpdateFrequency)
%
% In:
%   MatlabStream : name of the stream to create in the MATLAB environment (default: 'laststream')
%
%   DiskStream : DataRiver stream/file to read (default: '/tmp/DataRiver')
%
%   InputMetadata : Meta-data of the input stream. This is a struct or cell array of name-value
%                   pairs with meta-data fields to use. The mandatory fields are 'srate' and
%                   'chanlocs', where chanlocs is either a channel locations struct array, or
%                   a cell array of channel names, or the number of channels (in which case
%                   a cell array of the form {'A1','A2', ..., 'A32','B1','B2', ...} is created).
%                   Optionally, the field 'datasource' can be set to point to a dataset on disk or
%                   to a MATLAB workspace variable, to serve as the source of meta-data.
%
%   UpdateFrequency : update frequency, in Hz (default: 25)
%
% Examples:
%   % open a DataRiver input stream that is updated at 30 Hz and reads from the stream named 'C:\tmp\eeg'
%   run_readdatariver('UpdateFrequency',30,'DiskStream','C:\tmp\eeg','InputMetadata',{'srate',512,'chanlocs','eegchannels.sfp'});
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-11-19

% declare the name of this component (shown in the menu)
declare_properties('name','DataRiver stream');

default_prefix = fastif(ispc,'C:\tmp\','/tmp/');

% read options
arg_define(varargin, ...
    arg({'new_stream','MatlabStream'}, 'laststream',[],'New Stream to create. This is the name of the stream within the MATLAB environment.'), ...
    arg({'in_stream','DiskStream'}, [default_prefix 'DataRiver'],[],'Input DataRiver stream. This is the stream that shall be analyzed and processed.'), ...
    arg_sub({'in_metadata','InputMetadata'},{},@utl_parse_metadata, 'Meta-data of the input stream. These are fields as they appear in EEGLAB data sets; only sampling rate and channel labels are mandatory.'), ...
    arg({'update_freq','UpdateFrequency'},25,[],'Update frequency. New data is polled at this rate, in Hz.'));

% parse the meta-data spec
meta = utl_parse_metadata(in_metadata);

% load DataSuite, if necessary
global ds_lib;
if isempty(ds_lib)
    disp('Loading DataSuite...');
    startup_ds;
end

% open a DataRiver input stream (and make sure it gets cleaned up eventually)
datariver_stream = ds_OpenRead(in_stream);
meta.datariver_cleanup = onCleanup(@()ds_CloseRead(datariver_stream));

% create online stream stream (with all metadata from the dataset)
onl_newstream(new_stream,meta);

% create background reading job
onl_read_background(new_stream,@()read_block(datariver_stream),update_freq)


% background data reading function
function block = read_block(datariver_stream)
buffer = {};
while 1
    % try to get a new sample
    [received,sample] = ds_Read(datariver_stream);
    if ~received
        break; end
    % append it to the buffer
    buffer{end+1} = double(sample.Data(2:sample.nItems)');
end
block = [buffer{:}];

