function par_worker(port,portrange,timeout_heartbeat,varargin)
% Run a worker process on a cluster node.
% par_worker(Port,PortRange,TimeoutHeartbeat,Options...)
%
% Receives commands (string expressions) from the network, evaluates them, and sends off the result to 
% some collector (again as a string). Processing is done in a single thread. Matlab processes which
% are running this function can communicate with the function par_schedule, which schedules tasks
% to these workers.
%
% In:
%   Port: port number on which to listen for requests (default: 23547)
%         if the port is already in use, the next free one will be chosen, until the permitted
%         portrange is exceeded. if specified as 0, a free port is chosen directly
%
%         If you intend to share your worker with other users on the same network, it is a good idea
%         to keep the default port and portrange. Note that par_getworkers, however, will generally
%         determine the ports to use according to its own logic.
%
%   PortRange: number of successive ports to try, including the default/supplied port; (default: 0) 
%             * if 1, only the desired port will be tried
%             * if 0, the portrange will match the # of cores on the machine divided by the # of 
%               threads per MATLAB instance (parameter 'matlabthreads')
%
%   TimeoutHeartbeat : timeout for heartbeat messages; if nonzero and no such message is received 
%                      for the given period of time, the worker will terminate itself (default: 0)
%
%   Options... : optional name-value pairs, with possible names:
%                'backlog': backlog of queued incoming connections (default: 1)
%
%                'timeout_accept': timeout for accepting connections, in seconds (default: 3)
%
%                'timeout_recv': timeout for receiving data, in seconds (default: 5)
%
%                'timeout_send': timeout for sending results, in seconds (default: 10)
%
%                'timeout_dialout': timeout for dial-out connections (upon returning data), in seconds 
%                                   (default: 10)
%
%                'min_keepalive' : the minimum time for which the worker is kept alive after
%                                  startup, even if no heatbeat message is received at first, in
%                                  seconds (default: 300)
%
%                'retries_send': number of retries for sending a given result (default: 2)
%
%                'retry_wait': waiting period before re-trying, in seconds (default: 2)
%
%                'linger_send': time for which the connection is held open after completing a 
%                               send-back operation, in seconds (default: 3)
%
%                'update_check': information for code update checking (default: [] = no checking)
%                                this is a cell array of arguments to par_checkupdate
%
%                'matlabthreads': number of threads used by MATLAB on a machine (also determines
%                                 the portrange to use) (default: 4)
%
%                'token': the token used to indicate a heartbeat message (default: some random
%                         string)
%
% Notes:
%  * use multiple workers to make use of multiple cores
%  * use only ports that are not accessible from the internet
%  * request format: <task_id><collectoraddress_length><collectoraddress><body_length><body>
%    <task_id>: identifier of the task (needs to be forwarded, with the result,
%               to a some data collector upon task completion) (int)
%    <collectoraddress_length>: length, in bytes, of the data collector's address (int)
%    <collectoraddress>: where to send the result for collection (string, formatted as host:port)
%    <body_length>: length, in bytes, of the message body (int)
%    <body>: a MATLAB command that yields, when evaluated, some result in ans
%            (string, as MATLAB expression)
%            if an exception occurs, the exception struct (as from lasterror) replaces ans
%  * response format: <task_id><body_length><body>
%
% Examples:
%   % start a worker process on one of the default ports (and use a successively higher port number)
%   % if it is occupied, but try only as many ports as the computer has cores
%   par_worker;
%
%   % start a worker process on the default port, or return if that port is occupied
%   par_worker([],1);
%
%   % start a worker listening on some specific port, or return if that port is occupied
%   par_worker(51123,1);
%
%   % start a worker listening on some specific port, or one of the 5 successively higher ones
%   par_worker(51123,5);
%
%   % start a worker listening on some specific port, or one of the #cores successively higher ones
%   par_worker(51123);
%
%   % start a worker using default port settings, but customize some of the timing options
%   par_worker([],[],'timeout_accept',5,'timeout_recv',3,'timeout_send',15);
%
% See also:
%   par_schedule, par_beginschedule, par_endschedule, par_getworkers_ssh, par_getworkers_qsub
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-08-26

% Copyright (C) Christian Kothe, SCCN, 2010, christian@sccn.ucsd.edu
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

% this may be triggered during serialization quite frequently
warning off MATLAB:structOnObject

% parse arguments
if ~exist('port','var') || isempty(port)
    port = 23547; end
if ~exist('portrange','var') || isempty(portrange)
    portrange = 0; end
if ~exist('timeout_heartbeat','var') || isempty(timeout_heartbeat)
    timeout_heartbeat = 0; end
if ~isnumeric(port) || ~isscalar(port) || port-round(port) > eps
    error('The Port must be given as an integer.'); end
if ~isnumeric(portrange) || ~isscalar(portrange) || portrange-round(portrange) > eps
    error('The PortRange must be given as an integer.'); end

% read additional options
opts = hlp_varargin2struct(varargin, 'backlog',1, 'timeout_accept',3, 'timeout_dialout',10, 'timeout_send',10, ...
    'timeout_recv',5, 'min_keepalive',300, 'receive_buffer',64000, 'retries_send',2,'retry_wait',2, 'linger_send',3, 'update_check',[],'matlabthreads',4,'token','dsfjk45djf');

% make sure that the necessary Java/MATLAB code is present
if ~exist('ChunkReader','class')
    disp(['Trying to add Java class files in ' pwd]);
    javaaddpath(pwd);
end
if ~exist('ChunkReader','class')
    error('Cannot find the ChunkReader.class; please make sure that you have bundled it with the worker.'); end

% set portrange
if portrange == 0
    runtime = java.lang.Runtime.getRuntime();
    portrange = ceil(runtime.availableProcessors()/opts.matlabthreads);
end

try
    % try to adhere to the matlabthreads setting
    warning off MATLAB:maxNumCompThreads:Deprecated
    maxNumCompThreads(opts.matlabthreads);
catch
    disp('Note: could not restrict the # of threads to be used by MATLAB.');
end

% open a new server socket (trying the specified portrange)
opened = false;
for p = port : port+portrange-1
    try
        serv = ServerSocket(p, opts.backlog);
        opened = true;
        break;
    catch,end
end

% display status, or exit, if no port found
if opened
    disp(['This is ' char(InetAddress.getLocalHost()) '.']);
    disp(['Listening on port ' num2str(serv.getLocalPort()) '.']);
else
    disp('No free port found; exiting.');
    return; 
end

% set socket properties (e.g., making the function interruptible)
serv.setReceiveBufferSize(opts.receive_buffer);
serv.setSoTimeout(round(1000*opts.timeout_accept));
% make sure that the server socket will be closed when this function is terminated
cleaner = onCleanup(@()serv.close());

% try to open a keepalive listener socket, if desired
if timeout_heartbeat
    % this value is the time when we will terminate if no heartbeat is received before that time
    terminate_at = toc(uint64(0)) + max(opts.min_keepalive,timeout_heartbeat);
    try
        fprintf('Setting up heartbeat server... ');
        % we listen on the same port as the corresponding TCP port
        heartbeat = DatagramSocket(p);
        heartbeat.setSoTimeout(50);    % this is just the receive timeout (in ms)
        closer = onCleanup(@()heartbeat.close());
        disp('success.');
    catch e
        disp(['could not open heartbeat socket, exiting now. Reason: ' e.message]);
        return;
    end
end

tasknum = 1;
disp('waiting for connections...');
while 1
    try
        % wait for an incoming request
        conn = serv.accept();
        conn.setSoTimeout(round(1000*opts.timeout_recv));
        conn.setTcpNoDelay(1);
        disp('connected.');

        try
            % parse request
            in = DataInputStream(conn.getInputStream());
            cr = ChunkReader(in);
            taskid = in.readInt();
            collector = char(cr.readFully(in.readInt())');
            task = char(cr.readFully(in.readInt())');
            if ~isempty(opts.update_check)
                fprintf('received data (%.0fkb); checking for update first...\n',length(task)/1024);
                if par_haveupdate(opts.update_check{:})
                    disp('update available; terminating.');
                    conn.close();
                    return;
                else                    
                    disp('no update; now hand-shaking.');
                end
            else
                fprintf('received data (%.0fkb); hand-shaking.\n',length(task)/1024);
            end
            out = DataOutputStream(conn.getOutputStream());
            out.writeInt(taskid+length(collector)+length(task));
            out.flush();
            conn.close();
            
            % process task
            disp(['processing task ' num2str(taskid) ' (' num2str(tasknum) ') ...']); tasknum = tasknum+1;
            result = base64encode(par_evaluate(base64decode(task)));
            disp('done with task; opening back link...');
            
            for retry = 1:opts.retries_send
                try
                    % send off the result
                    idx = find(collector==':',1);
                    outconn = Socket();
                    destination = InetSocketAddress(collector(1:idx-1), str2num(collector(idx+1:end)));
                    outconn.connect(destination,round(1000*opts.timeout_dialout));
                    outconn.setTcpNoDelay(1);
                    outconn.setSoTimeout(round(1000*opts.timeout_send));
                    outconn.setSoLinger(true,3);
                    fprintf('connected; now sending (%.0fkb)...',length(result)/1024);
                    out = DataOutputStream(outconn.getOutputStream());
                    out.writeInt(taskid);
                    out.writeInt(length(result));
                    out.writeBytes(result);
                    out.flush();
                    fprintf('done; now closing...');
                    outconn.close();
                    disp('done.');
                    disp('waiting for connections...');
                    break;
                catch e
                    if isempty(strfind(e.message,'timed out'))
                        disp(['Exception during result forwarding: ' e.message]); end
                    disp('waiting before retry...');
                    pause(opts.retry_wait);
                    disp('retrying to open back link...');
                end
            end
            
        catch e
            conn.close();
            if ~isempty(strfind(e.message,'EOFException'))
                disp('cancelled by scheduler.');
            elseif isempty(strfind(e.message,'timed out'))
                disp(['Exception during task receive: ' e.message]); 
            end
        end
               
    catch e
        % accept timed out 
        if timeout_heartbeat
            % we're dependent on heartbeat messages; check what we got
            try
                while 1
                    % capture all datagram packets that have been enqueued so far
                    tmp = DatagramPacket(uint8(zeros(1,1024)),1024);
                    heartbeat.receive(tmp);
                    content = char(tmp.getData');
                    if strncmp(content,opts.token,length(opts.token))
                        % advance the terminate_at time
                        terminate_at = max(toc(uint64(0)) + timeout_heartbeat, terminate_at);
                    end
                end
            catch e
                % got an exception (= no more packets)
            end
            % check if the heartbeat timeout has expired
            if toc(uint64(0)) > terminate_at
                disp('No keep-alive message has been received in time; now terminating.');
                return;
            end
        end
        
        % check if we need to display an error message
        if isempty(strfind(e.message,'timed out'))
            disp(['Exception during accept: ' e.message]); end
    end
end
