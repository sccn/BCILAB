function result = par_fileclient(varargin)
% Allows to get files cached in a distributed in-memory filesystem.

args = arg_define(0:2,varargin, ...
    arg({'filename','FileName'},'',[],'Name of the file to operate on.'), ...
    arg({'operation','Operation'}, 'getfile', {'getfile','evict','list_hosts','list_files','terminate'}, 'Operation to perform.'), ...
    arg({'request_port','RequestPort'}, 34575, uint32([1 65535]), 'UDP port on which requests are being made.'), ...
    arg({'receive_buffer','ReceiveBuffer',},64000,[],'Size of the receive buffer.','guru',true), ...
    arg({'timeout_accept','TimeoutAccept',},300,[],'Timeout to start file transmission. In miliseconds. If a client takes longer than this to initiate transmission we give up.','guru',true), ...
    arg({'timeout_receive','TimeoutReceive',},60000,[],'Timeout to perform file transmission. In miliseconds. If a client takes longer than this to initiate transmission we give up.','guru',true), ...
    arg({'force','Force'}, false, [], 'Force process even if the file is not in cache.'), ...
    arg({'local_fallback','LocaLFallback'}, false, [], 'Fall back to local processing if unsuccessful.'), ...
    arg_nogui({'timeout_accept_granularity','TimeoutAcceptGranularity',},1000,[],'Timeout to start file transmission. In miliseconds. If a client takes longer than this to initiate transmission we give up.','guru',true), ...
    arg({'verbosity','Verbosity'}, 2, [], 'Verbosity level.'));

verbosity = args.verbosity;

% transmission timeout when force is set (i.e., we expect the operation to complete); in ms.
force_timeout = 15000;
% maximum amuont of memory we want to have in the java VM
max_java_memory = 2^26; % 64 MB

import java.io.*
import java.net.*
import java.lang.*

% create a request socket
if verbosity > 1
    fprintf('creating broadcast socket...'); end
sock = DatagramSocket();
%closer = onCleanup(@()sock.close());
if verbosity > 1
    fprintf('done.\n'); end

switch args.operation
    case 'getfile' % get a file from cache
        % local sanity checks
        if isempty(args.filename)
            error('Filename must be non-empty.'); end
        if ~exist(args.filename,'file')
            error('The given file "%s" does not exist.',args.filename); end
        % open a TCP server socket to be able to receive the result
        if verbosity > 1
            fprintf('creating TCP server socket...'); end
        serv = ServerSocket(0,1);
        serv.setReceiveBufferSize(args.receive_buffer);
        serv.setSoTimeout(quickif(args.force,force_timeout,round(args.timeout_accept_granularity)));
        serv_port = serv.getLocalPort();
        serv_addr = char(serv.getInetAddress().toString());
        closer = onCleanup(@()serv.close());
        if verbosity > 1
            fprintf('done.\nlistening on %s:%i.\n',serv_addr,serv_port); end
        if verbosity > 1
            fprintf('done.\n'); end
        % send the request packet
        verbose_send(verbosity,sock,sprintf('getfile\nfilename=%s\nforce=%i\nreturn_port=%i',args.filename,args.force,serv_port));
        % wait for the response
        try
            if verbosity > 1
                fprintf('waiting for incoming data...'); end
            t0 = tic;
            while toc(t0) < args.timeout_accept
                % note: this loop is mainly to allow us to interrupt it
                try
                    conn = serv.accept();
                    break;
                catch e
                    % accepted timed out; try again
                end
            end
            if verbosity > 1
                fprintf('got a connection.\nopening input stream...'); end
            % got one
            closer = onCleanup(@()conn.close());
            conn.setSoTimeout(quickif(args.force,force_timeout,round(args.timeout_receive)));
            % read the data
            in = DataInputStream(conn.getInputStream());
            if verbosity > 1
                fprintf('done.\nreceiving data'); end
            t0 = tic;
            full_length = round(in.readDouble());
            len = full_length;
            cr = ChunkReader(in);
            result = {};
            while len > max_java_memory
                result{end+1} = cr.readFully(max_java_memory);
                len = len - max_java_memory;
                if verbosity > 1
                    fprintf('.'); end
            end
            if len > 0
                result{end+1} = cr.readFully(len); end
            result = typecast(vertcat(result{:}),'uint8');
            if verbosity > 1
                fprintf('done (%.1f seconds; %.2fMB/s bandwidth).\n',toc(t0),length(result)/2^20/toc(t0)); end
        catch e
            % check if we need to display an error message
            if verbosity >= -1 && isempty(strfind(e.message,'timed out'))
                fprintf('exception during accept: %s\n',e.message);
            elseif verbosity > 0
                fprintf('connection timed out. No data obtained.');
            end
            if args.local_fallback
                fprintf('falling back to local file operation...');
                t0 = tic;
                f = fopen(args.filename,'r');
                filecloser = onCleanup(@()fclose(f));
                result = fread(f,inf,'*uint8');
                fclose(f);
                fprintf('done (%.1f seconds; %.2fMB/s bandwidth).\n',toc(t0),length(result)/2^20/toc(t0));
            else
                error('No result available.');
            end
        end
    case 'evict' % evict a file from cache
        verbose_send(verbosity,sock,sprintf('evict\nfilename=%s\nforce=%i',args.filename,args.force));
    case 'terminate' % terminate file hosts
        verbose_send(2,sock,sprintf('terminate\nforce=%i',args.force));
    case 'list_hosts'
        % TODO...
    otherwise
        error('Unsupported operation: %s',hlp_tostring(args.operation));
end

    % === utility functions ===
    
    function verbose_send(verbosity,sock,content)
        import java.net.*
        if verbosity > 1
            fprintf('sending request:\n%s\n',content); end
        sock.send(DatagramPacket(uint8(content),length(content),InetAddress.getByName('255.255.255.255'),args.request_port));
        if verbosity > 1
            fprintf('\n'); end
    end

end


