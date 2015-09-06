function varargout = par_testworker(address,func,varargin)
% Test a remote worker process by forwarding a task.
% Result = par_testworker(Address,Function,Arguments...)
% 
% In:
%   Address : network address of the worker (Host:Port)
%
%   Function : handle of the function to evaluate
%
%   Arguments... : optional arguments to the function
%
% Out:
%   Result : output of the function
% 
% Examples:
% 	par_testworker('192.168.0.1:23457',@sin,1000)
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-07-18

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

import java.io.*
import java.net.*
import java.lang.*


% some parameters
taskid = round(rand*100000);
disp(['assigning task id ' num2str(taskid)]);
receive_buffer_size = 64000;
receive_backlog = 10;
timeout_connect_ms = 3000;
timeout_send_ms = 5000;
timeout_recv_ms = 5000;
timeout_accept_ms = 5000;
use_hlpscope = true;
use_hlpgetresult = false;

if ischar(func)
    func = str2func(func); end;

% transfer the current scope (optional)
if use_hlpscope
    scope = hlp_resolveall;
    varargin = [{scope,func} varargin];
    func = @hlp_scope;
end

% construct message
fprintf('serializing task data...');
messagebody = hlp_serialize([{1},{func},varargin]);
disp('done.');


% --- set up TCP server to receive the reply ---

fprintf('opening TCP server to receive replies...');
while true
    try
        serv = ServerSocket(2000+round(rand*60000),receive_backlog);
        break;
    catch,end
end
serv.setReceiveBufferSize(receive_buffer_size);
serv.setSoTimeout(timeout_accept_ms);
close_server = onCleanup(@()serv.close());

% determine server's address for replies
return_port = num2str(serv.getLocalPort());
return_host = hlp_hostip;
return_address = [return_host ':' num2str(return_port)];
disp(['done; will receive on ' return_address]);


% --- send task ---

% determine target address
idx = find(address==':',1);
host = address(1:idx-1);
port = str2num(address(idx+1:end));
destination = InetSocketAddress(host,port);

% connect
fprintf('connecting to %s:%i...',host,port);
outconn = Socket();
outconn.connect(destination,timeout_connect_ms); % 3-second timeout
outconn.setTcpNoDelay(1);
outconn.setSoTimeout(timeout_send_ms);
outconn.setSoLinger(true,3);
close_outconn = onCleanup(@()outconn.close());
disp('done.');

% send message
fprintf('now sending (%.0fkb)...',length(messagebody)/1024);
out = DataOutputStream(outconn.getOutputStream());
out.writeInt(taskid);
out.writeInt(length(return_address));
out.writeBytes(return_address);
out.writeInt(length(messagebody));
out.write(messagebody);
out.flush();

% wait for handshake
fprintf('done; now waiting for checkum...');
in = DataInputStream(outconn.getInputStream());
checksum = in.readInt();
fprintf('done; checking...');
if checksum == taskid+length(return_address)+length(messagebody)
    fprintf('all clear');
else
    fprintf('ERROR: unexpected checksum');
end
% close socket
fprintf('; now closing...');
outconn.close();
disp('done.');


% --- receive reply ---

% wait for and get an incoming connection
fprintf('now waiting for reply...');
while true
    try
        conn = serv.accept();
    catch
        continue;
    end
    conn.setSoTimeout(timeout_recv_ms);
    conn.setTcpNoDelay(1);
    
    % receive reply
    fprintf('done; receiving data.\n');
    in = DataInputStream(conn.getInputStream());
    cr = ChunkReader(in);
    received_taskid = in.readInt();
    if received_taskid ~= taskid
        fprintf('got unexpected task id (%i instead of %i); back to waiting...\n',received_taskid,taskid);
        conn.close();
    else
        disp(['got correct task id: ' num2str(received_taskid)]);
    end
    fprintf('receiving reply...');
    replybody = typecast(cr.readFully(in.readInt()),'uint8');
    fprintf('got reply (%.0fkb).\n',length(replybody)/1024);
    
    fprintf('now closing...');
    conn.close();
    disp('done.');
    
    % construct result
    fprintf('parsing result...');
    reply = hlp_deserialize(replybody);
    disp('done.');
    
    % post-process result
    if isfield(reply{2},{'message','identifier','stack'})
        disp('Got an exception; traceback: ');
        rethrow(reply{2});
    else
        result = reply{2};
    end    
    break;
end

if use_hlpgetresult
    varargout = result;
else
    varargout = {result};
end