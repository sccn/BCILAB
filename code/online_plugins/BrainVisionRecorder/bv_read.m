% BV_READ  Reads BrainVision data over existing TCP socket.
%     DATA = BV_READ(H)
%     In:
%     h : handle to a BrainVision connection (previously opened via bv_open)
%     Out:
%     data : a block of new data, [#Channels x #Samples]

% Author: Hal Greenwald, The MITRE Corporation, 29-NOV-2011
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

function data = bv_read(h)
data = [];
if (~h.initialized)
    return;
end

header_size = 24;
try
    % check for existing data in socket buffer
    tryheader = pnet(h.handle, 'read', header_size, 'byte', 'network', 'view', 'noblock');
    while ~isempty(tryheader)
        
        % Read header of RDA message
        hdr = ReadHeader(h.handle);
        
        % Perform some action depending on the type of the data package
        switch hdr.type
            case 1       % Start message
                % Reset block counter to check overflows
                h.lastBlock = -1;
                pnet(h.handle, 'read', hdr.size - header_size);
                
            case 4       % 32Bit Data block
                % Read data and markers from message
                [datahdr, eegdata, markers] = ReadDataMessage(h.handle, hdr, h.channelCount);
                data = cat(2, data, eegdata);
                
                % check tcpip buffer overflow
                if h.lastBlock ~= -1 && datahdr.block > h.lastBlock + 1
                    disp(['******* Overflow with ' int2str(datahdr.block - h.lastBlock) ' blocks ******']);
                end
                h.lastBlock = datahdr.block;
                
                % Edit the following lines to process marker information
%                 for m = 1:datahdr.markerCount
%                     disp(markers(m));
%                 end
                
            case 3       % Stop message
                disp('Received STOP signal from Recorder');
                pnet(h.handle, 'read', hdr.size - header_size);
                bv_close(h);
                h.initialized = false;
                data = [];
                                                
            otherwise    % ignore all unknown types, but read the package from buffer
                pnet(h.handle, 'read', hdr.size - header_size);
        end
        if h.initialized
            tryheader = pnet(h.handle, 'read', header_size, 'byte', 'network', 'view', 'noblock');
        else
            tryheader = [];
        end
    end
catch er
    disp(er.message);
end

if ~isempty(data)
    % scale the data to uV (consistent with the pop_loadbv plugin)
    data = bsxfun(@times,data,h.resolutions');
end

%% ***********************************************************************
% Read the message header
function hdr = ReadHeader(con)
% con    tcpip connection object

% define a struct for the header
hdr = struct('uid',[],'size',[],'type',[]);

% read id, size and type of the message
% swapbytes is important for correct byte order of MATLAB variables
hdr.uid = pnet(con,'read', 16);
hdr.size = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
hdr.type = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));


%% ***********************************************************************
% Read a data message
function [datahdr, data, markers] = ReadDataMessage(con, hdr, channelCount)
% con       tcpip connection object
% hdr       message header
% datahdr   data header with information on datalength and number of markers
% data      data as one dimensional arry
% markers   markers as array of marker structs

% Define data header struct and read data header
datahdr = struct('block',[],'points',[],'markerCount',[]);

datahdr.block = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
datahdr.points = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
datahdr.markerCount = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));

% Read data in float format
data = swapbytes(pnet(con,'read', channelCount * datahdr.points, 'single', 'network'));
data = reshape(data,channelCount,datahdr.points);

% Define markers struct and read markers
markers = struct('size',[],'position',[],'points',[],'channel',[],'type',[],'description',[]);
for m = 1:datahdr.markerCount
    marker = struct('size',[],'position',[],'points',[],'channel',[],'type',[],'description',[]);
    
    % Read integer information of markers
    marker.size = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.position = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.points = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.channel = swapbytes(pnet(con,'read', 1, 'int32', 'network'));
    
    % type and description of markers are zero-terminated char arrays
    % of unknown length
    c = pnet(con,'read', 1);
    while c ~= 0
        marker.type = [marker.type c];
        c = pnet(con,'read', 1);
    end
    
    c = pnet(con,'read', 1);
    while c ~= 0
        marker.description = [marker.description c];
        c = pnet(con,'read', 1);
    end
    
    % Add marker to array
    markers(m) = marker;
end
