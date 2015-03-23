function par_filehost(varargin)
% This function starts a file host on the network.
% The request format accepted by this host on a UDP broadcast port (default 34575) is (without <>'s)
% <command>\n
% <argument1>=<value1>
% <argument2>=<value2>
% <argument3>=<value3>
% ...


% the database
persistent db;      % struct mapping of cache-tag --> data-blob
persistent lru;     % cell array of recently-used cache tags (sorted by increasing age)
if ~iscell(lru)
    lru = {}; end

fprintf('parsing arguments...');
settings = arg_define(varargin, ...
    arg({'all_hosts','AllHosts'}, {'computing-0-0.local','computing-0-1.local','computing-0-2.local','computing-0-3.local','computing-0-4.local','computing-0-5.local','computing-0-6.local','computing-0-7.local','computing-0-8.local','computing-0-9.local','computing-0-10.local'}, [], 'Hostnames where filehosts are running. Used for load-balancing.'), ...
    arg({'rank_override','RankOverride'}, 0, uint32([0 1000]), 'Rank of this instance (1-based). Only used if the current hostname is not in the set AllHosts.'), ...
    arg({'max_capacity','MaxCapacity'}, 140, [1 10 256 100000], 'Maximum memory capacity. In Gigabytes.'), ...
    arg({'accept_pattern','AcceptPattern'}, '*.*', [], 'File name pattern that is accepted by this host.'), ...
    arg({'per_pattern_capacity','PerPatternCapacity'}, {}, [], 'Per-pattern capacity. This is a cell array of {pattern,capacity,pattern,capacity,...} where pattern are wildcard patterns.','type','expression'), ...
    arg({'evict_policy','EvictPolicy'}, 'lru', {'lru','none'}, 'Eviction policy. What entries will be evicted from memory if capacity is exceeded. If none, no more files will be cached after capacity has been exceeded.'), ...
    arg({'request_port','RequestPort'}, 34575, uint32([1 65535]), 'UDP port on which requests are being received.'), ...
    arg({'timeout_backconnect','TimeoutBackconnect',},1000,[],'Timeout for back-connection. In miliseconds. If it takes longer than this to connect back to the requester we give up.','guru',true), ...
    arg({'timeout_transmit','TimeoutTransmit',},10000,[],'Timeout for transmission. In miliseconds. If we get a hiccup during back-transmission that is longer than this we give up.','guru',true), ...
    arg({'accept_packet_size','AcceptPacketSize'},1400,uint32([1 508 1472 65507]),'Maximum packet size accepted by this host.','guru',true));
fprintf('done.\n');

% maximum amount of memory used inside the Java VM
max_java_memory = 2^26; % 64 MB

fprintf('Importing java packages...');
import java.io.*
import java.net.*
import java.lang.*
fprintf('done.\n');

myhost = hlp_hostname;
fprintf('This is par_filehost on %s.\n',myhost);

% determine which rank we have
matching = strcmp(settings.all_hosts,myhost);
if any(matching)
    myrank = find(matching);
    fprintf('Rank of this host is %i/%i (based on AllHosts).\n',myrank,length(matching));
    fprintf('Allhosts was: %s\n',hlp_tostring(settings.all_hosts));
elseif settings.rank_override
    myrank = settings.rank_override;
    fprintf('Rank of this host is %i (based on RankOverride).\n',myrank);
else
    fprintf('Error: this host is not in AllHosts and no RankOverride was given. Exiting...\n');
    return;
end

% initialize request socket
fprintf('opening request socket...');
req_socket = DatagramSocket(settings.request_port,InetAddress.getByName('0.0.0.0')); %req_socket = DatagramSocket(settings.request_port,InetAddress.getByName('255.255.255.255'));
req_socket.setBroadcast(true);
req_socket.setSoTimeout(500);
closer = onCleanup(@()req_socket.close());
fprintf('done.\n');

while true
    try
        % get next request
        fprintf('\n\nwaiting for next packet...');
        packet = DatagramPacket(uint8(zeros(1,settings.accept_packet_size)),settings.accept_packet_size);
        while true
            try
                req_socket.receive(packet);
                break;
            catch e
                if isempty(strfind(e.message,'timed out'))
                    rethrow(e); end
            end
        end
        content = char(packet.getData');
        srcport = packet.getPort();
        srcaddr = packet.getAddress();
        fprintf('got packet:\n%s\n\n',content);
        
        % handle the request
        [command,args] = parse_request(content);
        fprintf('handling command %s; %i arguments.\n',command,length(fieldnames(args)));
        reply = [];
        switch command
            % --- primary commands ---
            case 'getfile'
                get_file(args);
            case 'evict'
                evict_file(args);
                % --- control commands ---
            case 'terminate'
                fprintf('  received termination command; exiting...\n');
                send_reply_udp(args,'ReplyContent',[myhost ' terminated']);
                return;
                % ---- diagnostics ---
            case 'getrank'
                reply = hlp_tostring(myrank);
            case 'allhosts'
                reply = hlp_tostring(settings.all_hosts);
            case 'hostname'
                reply = hlp_hostname;
            case 'max_capacity'
                reply = hlp_tostring(max_capacity);
            otherwise
                fprintf('  unrecognized command.\n');
        end
        % handle misc replies
        if ~isempty(reply)
            send_reply_udp(args,'ReplyContent',reply); end
    catch e
        fprintf('error: %s',hlp_handleerror(e));
        pause(2);
    end
end


    % === command handlers ===
    
    function get_file(varargin)
        import java.net.*
        import java.io.*
        % retrieve a file either from memory cache or from disk (if necessary)
        opts = hlp_varargin2struct(varargin, ...
            {'filename','FileName'},'', ...
            {'force','Force'},false, ...
            {'return_address','ReturnAddress'},'', ...
            {'return_port','ReturnPort'},'');
        
        % check if we're responsible, etc.
        [tag,namehash,responsible_host] = db_tag_if_responsible(opts.filename);
        fprintf('Handling get_file request for file "%s" (hash: %s), responsible_host=%i (%s)...',opts.filename,namehash,responsible_host,settings.all_hosts{responsible_host});
        if isempty(tag)
            fprintf('not responsible; ignoring.\n');
        else
            % got a request that we need to handle
            result = [];
            
            % check if the file matches our accept pattern
            pattern = ['^',strrep(strrep(settings.accept_pattern,'?','.'),'*','.{0,}'),'$'];
            if isempty(regexp(opts.filename,pattern, 'once'))
                fprintf('file does not match accept pattern; ignoring.\n');
            else
                if isfield(db, tag)
                    % record in cache: return it
                    fprintf('found record in cache.\n');
                    result = db.(tag);
                    lru = [{tag} lru(~strcmp(lru,tag))];
                else
                    % record not in cache
                    if ~exist(opts.filename,'file')
                        fprintf('file does not exist.\n');
                    else
                        % check if it would fit in our DB
                        capacity = settings.max_capacity*2^30;
                        filesize = getfield(dir(opts.filename),'bytes');
                        dbsize = getfield(whos('db'),'bytes');
                        fprintf('\ndb parameter: capacity:%.3fGB, filesize:%.3fGB, dbsize:%.3fGB\n',capacity/2^30,filesize/2^30,dbsize/2^30);
                        % determine if we have to load it
                        if (strcmp(settings.evict_policy,'lru') && (filesize < capacity)) || (strcmp(settings.evict_policy,'none') && (dbsize+filesize < capacity)) || opts.force
                            t0=tic; fprintf('loading file...');
                            f = fopen(opts.filename,'r');
                            filecloser = onCleanup(@()fclose(f));
                            result = fread(f,inf,'*uint8');
                            fprintf('%.1f seconds (read bandwidth: %.2fMB/s).\n',toc(t0),filesize/2^20/toc(t0));
                            % cache it if applicable
                            if strcmp(settings.evict_policy,'none') && (dbsize+filesize < capacity)
                                fprintf('committing to db.\n');
                                db.(tag) = result;
                                lru = [{tag} lru(~strcmp(lru,tag))];
                            elseif strcmp(settings.evict_policy,'lru') && (filesize < capacity)
                                % evict least-recently used records until we're within the capacity limits
                                fprintf('committing to db.\n');
                                db.(tag) = result;
                                lru = [{tag} lru(~strcmp(lru,tag))];
                                if filesize+dbsize > capacity
                                    fprintf('evicting least-recently used cache records');
                                    initialsize = getfield(whos('db'),'bytes');
                                    while ~isempty(fieldnames(db))
                                        fprintf('.');
                                        db = rmfield(db,lru(end));
                                        lru(end) = [];
                                        if getfield(whos('db'),'bytes') < capacity
                                            break; end
                                    end
                                    fprintf('done. Evicted %.3fGB of data.\n',(initialsize - getfield(whos('db'),'bytes'))/2^30);
                                end
                            end
                        else
                            fprintf('not loading.\n');
                        end
                    end
                end
            end
            
            % send off the result
            if isa(result,'uint8')
                [return_address,return_port] = get_return_endpoint(opts);
                fprintf('preparing to return result to requestor at %s:%i\n',char(return_address.toString()),return_port);
                % send the data back
                fprintf('opening back-connection to transmit data...');
                outconn = Socket();                
                destination = InetSocketAddress(return_address, return_port);
                outconn.connect(destination, round(settings.timeout_backconnect));
                sockcloser = onCleanup(@()outconn.close());
                outconn.setSoTimeout(round(settings.timeout_transmit));
                outconn.setSoLinger(true,3);
                fprintf('connected.\nnow sending (%.0fkb)',length(result)/1024);
                t0 = tic;
                out = DataOutputStream(outconn.getOutputStream());
                out.writeDouble(length(result));
                numsplits = ceil(length(result)/max_java_memory);
                for i=0:numsplits-1
                    range = 1+floor(i*length(result)/numsplits) : min(length(result),floor((i+1)*length(result)/numsplits));
                    out.write(typecast(result(range),'int8'),0,length(range));
                    fprintf('.');
                end
                out.flush();                
                fprintf('done (%.1f seconds; %.2fMB/s bandwidth)\nnow closing...',toc(t0),length(result)/2^20/toc(t0));
                fprintf('done.\ntransaction complete.\n');
            end
        end
    end

    function evict_file(varargin)
        % evict a file record from cache
        opts = hlp_varargin2struct(varargin, ...
            {'filename','FileName'},'');
        [tag,namehash,responsible_host] = db_tag_if_responsible(opts.filename);
        fprintf('Handling evict request for file "%s" (hash: %s), responsible_host=%i (%s)...',opts.filename,namehash,responsible_host,settings.all_hosts{responsible_host});
        if isfield(db,tag)
            db = rmfield(db,tag);
            lru(strcmp(lru,tag)) = [];
            fprintf('evicted record.\n');
        elseif ~isempty(tag)
            fprintf('no record in database.\n');
        else
            fprintf('no record in database and not responsible.\n');
        end
    end

    % === helper functions ===
    
    function send_reply_udp(varargin)
        import java.net.*
        % send a reply message over UDP
        opts = hlp_varargin2struct(varargin, ...
            {'reply_content','ReplyContent'},'');
        [return_address,return_port] = get_return_endpoint(opts);
        fprintf('sending reply packet to %s:%i...',char(return_address.toString()),return_port);
        repl_sock = DatagramSocket();
        repl_packet = DatagramPacket(uint8(opts.reply_content),length(opts.reply_content),return_address,return_port);
        repl_sock.send(repl_packet)
        fprintf('done.\n')
    end

    function [return_address, return_port] = get_return_endpoint(varargin)
        import java.net.*
        % determine the return endpoint based on the global srcaddr/port or override options
        opts = hlp_varargin2struct(varargin, ...
            {'return_address','ReturnAddress'},'', ...
            {'return_port','ReturnPort'},'');
        % allow overriding the return address / port
        fprintf('determining return address/port...');
        if isfield(opts,'return_address') && ~isempty(opts.return_address)
            return_address = InetAddress.getByName(opts.return_address);
        else
            return_address = srcaddr;
        end
        if isfield(opts,'return_port') && ~isempty(opts.return_port)
            if ischar(opts.return_port)
                return_port = str2num(opts.return_port);
            else
                return_port = opts.return_port;
            end
        else
            return_port = srcport;
        end
        fprintf('done.\n');
    end

    function [tag,namehash,responsible_host] = db_tag_if_responsible(filename)
        % returns a database tag if we're responsible for the given filename or empty otherwise
        % get a hash of the filename
        namehash = hlp_cryptohash(filename);
        % check if we're responsible for this file
        responsible_host = 1+mod(hex2dec(namehash(1:5))+1,length(settings.all_hosts));
        if responsible_host == myrank
            tag = ['x_' namehash];
        else
            tag = [];
        end
        
    end

    function [command,args] = parse_request(content)
        content(content==0) = [];
        % parse a network request packet
        lines = hlp_split(content,sprintf('\n'));
        % first line is the command
        command = strtrim(lower(lines{1}));
        % followed by one argument per line
        args = [];
        for l=2:length(lines)
            line = strtrim(lines{l});
            % strip comments (beginning with ;)
            comment = find(line == ';',1);
            if ~isempty(comment)
                line = line(1:comment-1); end
            if isempty(line)
                continue; end
            % handle assignments
            equals = find(line == '=',1);
            if isempty(equals)
                fprintf('the line "%s" is not a valid argument declaration; ignoring.\n',line); end
            key = strtrim(line(1:equals-1));
            val = strtrim(line(equals+1:end));
            try
                args.(key) = val;
            catch
                fprintf('argument "%s" is not a valid field name; ignoring.\n',key);
            end
        end
    end
end

