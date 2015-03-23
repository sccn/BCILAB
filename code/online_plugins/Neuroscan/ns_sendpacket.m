function ns_sendpacket(h,id,code,request,bodysize)
% Prepare and send a data packet to Neuroscan Scan server
% ns_sendpacket(h,id,code,request,bodysize)
%     
%
% In:
%   h: PNET connection handle
%
%   id: Neuroscan message ID string
%
%   code: Neuroscan message code value
%
%   request: Neuroscan message request value
%
%   bodysize: Size of message body
%
% Note: See the Neuroscan Scan documentation for description of message ID,
%   code, and request values
%
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald


% assemble the packet
packet = zeros(12,1,'uint8');
packet(1:4) = id;
packet(5:6) = fliplr(typecast(int16(code),'int8'));
packet(7:8) = fliplr(typecast(int16(request),'int8'));
packet(9:12) = fliplr(typecast(int32(bodysize),'int8'));
% write the packet
pnet(h,'write', packet)
    
end
    