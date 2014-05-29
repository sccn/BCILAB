function ns_close(h)
% Close a TCP connection to Neuroscan Scan
% ns_close(h)
%
%
% In:
%   h : handle to an existing Neuroscan connection
%     
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald


% send message to Scan to stop sending data
ns_sendpacket(TCP_Connection,'CTRL',3,4,0);
% and indicate that connection is closing
ns_sendpacket(TCP_Connection,'CTRL',1,2,0);

if ~h.initialized
    return;
end
disp('Cleaning up connection to Neuroscan Scan');
pnet(h.handle, 'close');
if evalin('base', ['exist(' 'sprintf(''%s'', h.name)' ')'])
    evalin('base',['clear ' sprintf('%s',h.name) ';']);
end