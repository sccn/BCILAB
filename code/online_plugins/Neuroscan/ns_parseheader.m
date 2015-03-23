function hdr = ns_parseheader(header)
% Parse the header info packet returned by Neuroscan Scan server
% hdr = ns_parseinfo(dataBytes)
%
%
% In:
%   dataBytes : 12 byte array returned by Neuroscan server as header
%
% Out:
%   hdr: Structure containing each element of packet header as a field.
%        Header information includes:
%
%       id: ID string, 'CTRL', 'FILE', or 'DATA'
%
%       code: Control code, 1 (General), 2 (Server), or 3 (Client)
%
%       req: Request value, see the Neuroscan Acquire manual
%
%       bodysize: Size (in bytes) of message attached to header,
%           or 0 if packet does not contain a data body (for example start
%           or stop acquisition messages)
%
%     
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald


hdr = struct('id',[],'code',[],'req',[],'bodysize',[]);

hdr.id = char(header(1:4));
hdr.code = double(typecast(fliplr(uint8(header(5:6))), 'uint16'));   
hdr.req = double(typecast(fliplr(uint8(header(7:8))), 'uint16'));
hdr.bodysize = double(typecast(fliplr(uint8(header(9:12))), 'uint32'));

   