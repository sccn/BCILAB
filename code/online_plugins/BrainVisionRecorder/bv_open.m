% BV_OPEN  Opens a TCP connection to BrainVision Recorder
%     H = BV_OPEN
%     Out:
%     h : handle to a newly opened BrainVision connection
%     Connection is automatically closed when handle is deleted

% Author: Hal Greenwald, The MITRE Corporation, 28-OCT-2011
% This software incorporates portions of the BrainVision Recorder RDA Client.  Used with permission of Brain Products GmbH.
   
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

function h = bv_open(hostname)

% Establish connection to BrainVision Recorder Software 32Bit RDA-Port
% (use 51234 to connect with 16Bit Port)
h.handle = pnet('tcpconnect', hostname, 51244);
h.initialized = false;

% Check established connection and display a message
stat = pnet(h.handle,'status');
if stat > 0
    disp('BrainVision Recorder connection established');
else
    disp('BrainVision Recorder connection failed');
    return;
end

%Wait for start message
disp('Waiting for start message from Recorder...');
header_size = 24;
while (~h.initialized)
    try
        data = pnet(h.handle, 'read', header_size, 'byte', 'network', 'view', 'noblock');
        if (~isempty(data))
            hdr = ReadHeader(h.handle);
            if (hdr.type == 1)
                disp('Starting online data collection...');
                [h.channelCount, h.samplingInterval, h.resolutions, h.channelNames] = ReadStartMessage(h.handle, hdr);
                h.initialized = true;
                h.lastBlock = -1;
                h.cleanup = onCleanup(@()pnet(h.handle, 'close'));
            else
                pnet(h.handle, 'read', hdr.size - header_size);
            end
        end
    catch er
        disp(er.message);
        return;
    end
end

function hdr = ReadHeader(con)
% con    tcpip connection object

% define a struct for the header
hdr = struct('uid',[],'size',[],'type',[]);

% read id, size and type of the message
% swapbytes is important for correct byte order of MATLAB variables
% pnet behaves somehow strange with byte order option
hdr.uid = pnet(con,'read', 16);
hdr.size = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
hdr.type = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));

function [channelCount, samplingInterval, resolutions, channelNames]  = ReadStartMessage(con, hdr)
% con    tcpip connection object
% hdr    message header

% read EEG properties
channelCount = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
samplingInterval = swapbytes(pnet(con,'read', 1, 'double', 'network'));
resolutions = swapbytes(pnet(con,'read', channelCount, 'double', 'network'));
allChannelNames = pnet(con,'read', hdr.size - 36 - channelCount * 8);
channelNames = SplitChannelNames(allChannelNames);

function channelNames = SplitChannelNames(allChannelNames)
% allChannelNames   all channel names together in an array of char
% channelNames      channel names splitted in a cell array of strings

% cell array to return
channelNames = {};

% helper for actual name in loop
name = [];

% loop over all chars in array
for i = 1:length(allChannelNames)
    if allChannelNames(i) ~= 0
        % if not a terminating zero, add char to actual name
        name = [name allChannelNames(i)];
    else
        % add name to cell array and clear helper for reading next name
        channelNames = [channelNames {name}];
        name = [];
    end
end
