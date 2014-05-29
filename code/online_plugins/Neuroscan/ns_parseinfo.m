function basicinfo = ns_parseinfo(dataBytes)
% Parse the basic info packet returned by Neuroscan Scan
% basicinfo = ns_parseinfo(dataBytes)
%
%
% In:
%   dataBytes : 28 byte array returned by Neuroscan server after requesting
%   basic info
%
% Out:
%   basicinfo: Structure containing each element of basic EEG info as a
%   field. EEG info includes:
%
%       size: size of info array (in bytes)
%
%       numChan: number of EEG data channels
%
%       numEventChan: number of event marker channels
%
%       samplesPerBlock: the number of data samples transmitted each block
%
%       srate: data sampling rate in Hz
%
%       bytesPerSamples: number of bytes per sample, either 2 (16 bits
%       per sample) or 4 (32 bits per sample)
%
%       resolution: the value in microvolts represented by the least
%       significant bit 
%
%     
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald

% basic EEG info 
basicinfo = struct('size',[],'numChan',[],'numEventChan',[],'samplesPerBlock',[],...
    'srate', [], 'bytesPerSample', [], 'resolution', []);

basicinfo.size = double(typecast((uint8(dataBytes(1:4))), 'int32')); 
basicinfo.numChan = double(typecast((uint8(dataBytes(5:8))), 'int32')); 
basicinfo.numEventChan = double(typecast((uint8(dataBytes(9:12))), 'int32'));
basicinfo.samplesPerBlock = double(typecast((uint8(dataBytes(13:16))), 'int32'));
basicinfo.srate = double(typecast((uint8(dataBytes(17:20))), 'int32'));
basicinfo.bytesPerSample = double(typecast((uint8(dataBytes(21:24))), 'int32'));
basicinfo.resolution = double(typecast((uint8(dataBytes(25:28))), 'single'));

end
