% RUN_READBRAINVISION  Receives real-time data from BrainVision Recorder
%     RUN_READBRAINVISION(MATLABSTREAM, CHANNELRANGE, SAMPLINGRATE, UPDATEFREQUENCY)
%     In:
%     StreamName : name of the stream to create in the MATLAB environment (default: 'laststream')
%     ChannelRange : numeric vector of channel indices that should be recorded   
%     SamplingRate : sampling rate for the amplifier, in Hz (default: 256)
%     UpdateFrequency : update frequency, in Hz (default: 25)
%     
%     Examples:
%     % open an input stream that is sampled at 512 Hz and updated at 30 Hz
%     % (using the default name for the stream)
%     run_readbrainvision('UpdateFrequency',30, 'SamplingRate',512);

% Author: Hal Greenwald, The MITRE Corporation, 29-NOV-2011
   
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function run_readbrainvision(varargin)

persistent bvhostname;
if isempty(bvhostname)
    bvhostname = '127.0.0.1';
end
persistent bvstreamname;
if isempty(bvstreamname)
    bvstreamname = 'laststream';
end
persistent bvfreq;
if isempty(bvfreq)
    bvfreq = 10;
end

% declare the name of this component (shown in the menu)
declare_properties('name','BrainVision Recorder');

% read options
opts = arg_define(varargin, ...
    arg({'src_hostname','InputHost'}, bvhostname,[],'Source TCP hostname. Can be a computer name, URL, or IP address.'), ...
    arg({'new_stream','StreamName','MatlabStream'}, 'laststream',[],'MATLAB Stream Name. A variable with this name will be created in the MATLAB workspace to hold the stream''s data. If such a variable already exists it will be overridden.'), ...
    arg({'update_freq','UpdateFrequency'}, bvfreq,[0 Inf],'Update frequency. New data is polled at this rate, in Hz.'));

bvhostname = opts.src_hostname;
bvstreamname = opts.new_stream;
bvfreq = opts.update_freq;

% open a connection
conn = bv_open(opts.src_hostname);
if (~conn.initialized)
    return;
end
conn.name = opts.new_stream;

%Create and initialize online stream
onl_newstream(opts.new_stream, 'srate', 1e6/conn.samplingInterval, 'chanlocs', conn.channelNames, 'data', zeros(length(conn.channelNames),0,0),'xmin',toc(uint64(0)));

% start background acquisition
onl_read_background(opts.new_stream,@()bv_read(conn), opts.update_freq);