function par_broadcast(command)
% Broadcast a command to all registered workers.
%
% This function can be used to change the state of all workers while they are running; for instance,
% when editing parallel code it is sometimes useful to change a function and have all workers update
% their (cached) version of the function without restarting them, or without clearing other caches
% on the workers. This can be accomplished by running par_broadcast('clear myfunction')
%
% Note that, when a worker is busy it will not be reachable by this function and you will get an
% error.
%
% In:
%   Command : the command to execute
%
%                                Christian Kothe, Syntrogi
%                                2015-04-16

% Copyright (C) Christian Kothe, Syntrogi, 2015, christian.kothe@syntrogi.com
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

timeout_dialout = 3;
timeout_send = 2;

pool = par_globalsetting('pool');
if isempty(pool)
    fprintf('No workers connected.\n');
else
    fprintf('Issueing command...\n');
    for p=1:length(pool)
        fprintf('connecting to %s...',pool{p});
        atsign = pool{p}=='@';
        if any(atsign)
            pool{p} = pool{p}(find(atsign,1)+1:end); end
        colon = find(pool{p}==':',1);
        conn = Socket();
        destination = InetSocketAddress(pool{p}(1:colon-1), str2num(pool{p}(colon+1:end)));
        conn.connect(destination,round(1000*timeout_dialout));
        conn.setTcpNoDelay(1);
        conn.setSoTimeout(round(1000*timeout_send));
        conn.setSoLinger(true,3);
        fprintf('connected; waiting for ready-to-send...');
        in = DataInputStream(conn.getInputStream());
        if in.readInt()==12345
            fprintf('confirmed.\nnow sending...');
            out = DataOutputStream(conn.getOutputStream());
            % encode the task
            taskid = round(23+rand()*43535);
            task = {taskid,@eval,command}; %#ok<CCAT>
            task = fast_encode(hlp_serialize(task));
            out.writeInt(taskid);       % task id
            out.writeInt(0);            % length of return address (none)
            out.writeInt(length(task)); % length of task
            out.writeBytes(task);       % task description
            out.flush();
            if in.readInt()==taskid+length(task)
                fprintf('confirmed.\n');
            else
                fprintf('failed. The worker likely will not process the command.\n');
            end
        else
            fprintf('failed.\nskipping...\n');
        end
        fprintf('closing connection...');
        conn.close();
        fprintf('done.\n');
    end
end
