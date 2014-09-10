function h = ns_open(hostname, port)
% Open a TCP connection to Neuroscan Recorder
% h = ns_open(Hostname, Port)
%     
%
% In:
%   Hostname: Source TCP hostname. Can be a computer name, URL, or IP
%       address
%
%   Port : the port on which to connect to the TCP host
%
% Out:
%   h : handle to a newly opened Neuroscan connection
%
% Author: Visual Attention and Cognition Lab, Dan Roberts, and Nick Peñaranda, George Mason University, Spring 2014
%         Released under the GPLv3, see COPYING.txt
%         Based on the BrainVision BCILAB plug-in by Hal Greenwald


% open the connection
h.handle = pnet('tcpconnect', hostname, port);
ConnectionStatus = pnet(h.handle,'status');

if ConnectionStatus > 0
    disp('Neuroscan Scan connection established');
else
    error('Neuroscan Scan connection failed - check the IP and port');
end

% flush the buffer
bufferData = pnet(h.handle,'read', 'noblock');
while ~isempty(bufferData);
    bufferData = pnet(h.handle,'read', 'noblock');
end

% request basic header
ns_sendpacket(h.handle,'CTRL',3,5,0);

% try
    % read reader and basic info
    packetBytes = pnet(h.handle,'read', 40, 'uint8');
    
    % parse header data
    header = ns_parseheader(packetBytes(1:12));
    basicinfo = ns_parseinfo(packetBytes(13:end));

    % attach Neuroscan parameters to connection handle
    infoFields = fields(basicinfo);
    for f = 1:length(infoFields)
        h.(infoFields{f}) = basicinfo.(infoFields{f});
    end
    
    h.totalChan = h.numChan + h.numEventChan;
    
    if h.numEventChan ~= 0
        h.markerChanIdx = h.numChan + 1;
    end
    
    % number of bytes for neuroscan header
    h.headerSize = 12;
    
    % save the block size, equal to the number of channels (including
    % marker channel) x the number of samples per block x the number of
    % bytes per sample
    h.dataBlockSize = (basicinfo.numChan + basicinfo.numEventChan) ...
        * basicinfo.samplesPerBlock * basicinfo.bytesPerSample;
    
    % save the data type (16 or 32 bit integer) for later casting
    if h.bytesPerSample == 2
        h.datatype = 'int16';
    elseif h.bytesPerSample == 4
        h.datatype = 'int32';
    else
        error('expecting either 2 or 4 bytes per sample');
    end

    h.cleanup = onCleanup(@()pnet(h.handle, 'close'));
        
    % instruct Scan to begin sending data
    ns_sendpacket(h.handle,'CTRL',3,3,0);    
    h.initialized = true;

% catch er
%     disp(er.message);
%     return;
% end

end

