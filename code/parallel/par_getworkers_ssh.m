function pool = par_getworkers_ssh(varargin)
% Acquire workers on some remote machines and return hostnames and ports of those that are available.
% Pool = par_getworkers_ssh(Hostnames,ProcessorsPerNode,MatlabThreads,StartupCommand,ShareResources,SharingPort,LoggingPath,IdentityFile,BinaryWorker,BinaryNmae,MatlabCompilerRoot,MinFreeMemory,MaxCPULoad,LoadWindow,AvoidProcess,RecruitMachines,VerboseOutput)
%
% This function attempts to ensure that the desired number of worker processes on remote machines
% is available for participating in distributed computations. Such computations are usually being
% submitted in the form of tasks to the respective worker processes, via the function par_schedule.
%
% Alternatively, the workers can also be started manually on the remote machines (this is the 
% recommended way if Windows machines are involved or use of SSH is not permitted) and their names
% would then just be listed in the cluster settings (rather than to configure acquisition behavior).
%
% In:
%   --- what shall be acquired ---
%
%   Hostnames : cell array of hostnames on which to acquire workers (duplicates and port assignments
%               will be ignored).
%
%   ProcessorsPerNode : number of processes to acquire per node; 0 means as many as there are cores
%                       on each node, divided by the value of the 'matlabthreads' parameter
%                       (default: 0)
%
%   MatlabThreads : assumed number of threads per MATLAB instance (default: 4)
%
%   --- how shall it be acquired ---
%
%   StartupCommand : startup command, if any; if binary_worker is true, this is a sequence of shell 
%                    commands, otherwise it is a sequence of MATLAB statements which contain the
%                    placeholder %d where the TCP port # would be inserted; If the startup_command
%                    is passed as 'autogenerate', a bcilab-specific startup command will be
%                    generated (default: 'autogenerate')
%
%   ShareResources : whether worker that are already running will be returned as part of the worker 
%                    set, and are thus potentially shared with other users. If false, the requested
%                    number of workers will always be started and be private to the requester -
%                    however, potentially over-subscribing the node. (default: true)
%
%   SharingPort : the starting port for shared workers; if 0, this is a hash of the startup_command 
%                 (default: 0)
%
%   --- associated files ---
%
%   LoggingPath : default logfile path; make sure that this does not conflict with other users' log 
%                 paths. (default: '')
%
%   IdentityFile : SSH identity file, if needed (implying the -i option for ssh) (default: '')
%   
%   --- optional binary worker support ---
%
%   BinaryWorker : whether to use assume a binary worker implementation rather than the MATLAB 
%                  worker (default: false)
%
%   BinaryName : name of the worker binary (default: 'build')
%
%   MatlabCompilerRoot : installation path of the MATLAB compiler runtime; if empty, a few locations 
%                        will be searched for an MCR that corresponds to the MATLAB version that is
%                        running this function. If nothing is found, it is assumed that the MCR has
%                        been installed into the system path  and is found automatically. (default: '')
%
%   --- exclusion criteria ---
%
%   MinMemory : minimum memory (in bytes) that needs to be available for an instance to be started 
%               (default: 4GB)
%
%   MaxCPULoad : maximum acceptable average CPU load (summed over all cores) for a node to be 
%                considered (default: Inf)
%
%   LoadWindow : CPU load measurement window, in minutes; can be 1, 5 or 15 (default: 15)
%
%   AvoidProcess : avoid nodes with the given process running (default: 'avoidme')
%
%   ShutdownTimeout : if non-zero, the worker will shut itself down if it has not received a
%                     heartbeat signal from a client in the given time, in seconds. (default: 0)
%
%   RecruitMachines : whether to recruit other machines, rather than just list them (default: true)
%
%
%   --- misc ---
%
%   VerboseOutput : whether to show verbose outputs during acquisition (default: true)
%
%
% Out:
%   Pool : cell array of 'hostname:port' strings specifying the list of available machines
%          Note: in the case that no return value is requested, the global variable 
%                tracking.parallel.pool will receive this result. This is the recommended way to
%                use par_getworkers_ssh, as par_schedule uses this pool by default.
%
% Notes:
%   If workers are shared between users, code which maintains global state (in particular
%   user-specific configuration options) may show unexpected behavior, as these variables are not
%   duplicated across users.
%
%   If resources are not shared, it must be considered that each call to par_getworkers_ssh will launch
%   another set of worker processes on the target nodes; there is no mechanism in place to
%   automatically terminate the remote workers when the MATLAB session ends.
%
%   This function works only across POSIX-compliant (Unix/Linux/Mac) systems. On standard Windows
%   workstations, the workers must be started manually (unless POSIX emulation software is
%   installed).
%
% Examples:
%   % acquire workers on the three hosts (assuming bcilab is in the same file-system location remotely as here)
%   par_getworkers_ssh({'machine1','machine2','machine3'})
%
%   % try to start fast pre-compiled workers (assuming bcilab is in the same file-system location remotely as here)
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'}, 'binary_worker',true)
%
%   % MATLAB version again, also setting a custom logging path
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'}, 'logging_path','/home/christian/workerlogs/')
%
%   % like before, but this time only consider machines that have at least 8GB of free RAM and at most 10% CPU load
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'}, 'logging_path','/home/christian/workerlogs/', 'min_memory',2^33, 'max_cpuload',0.1)
%
%   % like before, but measuring the CPU load within the last minute (instead to 15 minutes)
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'}, 'logging_path','/home/christian/workerlogs/', 'min_memory',2^33, 'max_cpuload',0.1, 'load_window',1)
%
%   % acquire workers on the three hosts, this time using a custom startup command (the %d is where the port goes)
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'},'startup_command','cd path/to/bcilab/; bcilab(''myconfig'',''worker'',{%d,1}')
%
%   % try to start fast pre-compiled workers (instead of entire MATLAB sessions); note that the 
%   % startup command is this time a UNIX shell command which CD's into the directory where the run_build.sh is located
%   par_getworkers_ssh('hostnames', {'machine1','machine2','machine3'}, 'binary_worker',true, 'startup_command','cd path/to/bcilab/build/distrib')
%
% See also:
%   par_worker
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-02-15

% Copyright (C) Christian Kothe, SCCN, 2011, christian@sccn.ucsd.edu
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

global tracking;

arg_define(varargin,...
    arg({'hostnames','Hostnames'},{''},[],'Host names to acquire. Host names or IP addresses of machines on which to consider starting workers. If this is left empty, the current global worker pool will be used.'), ...
    arg({'processors_per_node','ProcessorsPerNode'},0,[],'Processors per node. A value of 0 means that as many processors as there are cores on each node shall be used, but divided by the value of the "matlabthreads" (Number of MATLAB threads) parameter.'), ...
    arg({'matlabthreads','MatlabThreads'},4,[],'Number of MATLAB threads. This is the number of threads that each worker uses internally.'), ...
    arg({'startup_prefix','StartupPrefix'},'',[], 'Startup prefix. Any startup lines to run before running the main BCILAB command.'),...
    arg({'startup_command','StartupCommand'},'autogenerate',[], 'Startup command. If the binary_worker parameter (Use compiled workers) is true, this is a sequence of shell commands, otherwise it is a sequence of MATLAB statements which contain the placeholder %d where the TCP port # would be inserted; If set to ''autogenerate'', a startup command specific to the current BCILAB environment will be generated.'),...
    arg({'share_resources','ShareResources'},true,[],'Share resources. Whether workers that are already running (started by other users) will be returned as part of the worker set, and are thus potentially shared with other users. If false, the requested number of workers will always be started and be private to the requester -- however, potentially over-subscribing the capacity of the node.'),...
    arg({'sharing_port','SharingPort'},0, [], 'Port for sharing. The lowest port used by shared workers -- if 0, this is computed based on the startup_command (recommended).'), ...
    arg({'logging_path','LoggingPath'},env_translatepath('home:/.bcilab/logs/workers'), [], 'Logging path. Make sure that this does not conflict with other users'' log paths.'), ...
    arg({'identity_file','IdentityFile'},'', [], 'SSH identity file. Optional, if needed for passwordless login. This corresponds to the -i option in ssh.'), ...
    arg({'binary_worker','BinaryWorker'},false,[], 'Use compiled workers. Whether to start a given binary worker implementation rather than the MATLAB worker.'), ...
    arg({'binary_name','BinaryName'},'build', [], 'Name of the binary. The name of the BCILAB binary (if using compiled workers).'), ...
    arg({'matlab_command','MATLABCommand'},'matlab', [], 'Command to start MATLAB. This is the command to run on the command line to start the desired version of MATLAB.'), ...
    arg({'mcr_root','MatlabCompilerRoot'},'', [], 'MATLAB Compiler directory. Installation path of the MATLAB compiler runtime; if empty, a few locations will be searched for an MCR that corresponds to the MATLAB version that is running this function (this assumes that the compiler is installed in the same directory remotely as it is locally). If nothing is found, it is assumed that the MCR has been installed into the system path and is found automatically.'), ...
    arg({'min_memory','MinFreeMemory'},2^32,[], 'Minimum free memory. Minimum amount of RAM memory (in bytes) that needs to be available for an instance to be started.'),...
    arg({'max_cpuload','MaxCPULoad'},Inf, [], 'Maximum CPU load. The maximum acceptable average CPU load (summed over all cores) for a node to be considered.'), ...
    arg({'load_window','LoadWindow'},15,[],'Load estimation window. The time window in seconds (into the past) over which the CPU load will be estimated to determine whether a worker shall be started on a given node.'), ...
    arg({'avoid_proc','AvoidProcess'},'avoidme',[],'Process to avoid. Optionally a process name that indicates that no workers shall be spawned if the process is running on some machine (e.g. a time-critical computation).'),...
    arg({'shutdown_timeout','ShutdownTimeout'},0,[],'Shutdown timeout. If non-zero, the worker will be shut down if it has not received a heartbeat signal from a client in the given time frame.'),...
    arg({'recruit_machines','RecruitMachines'},true,[],'Recruit workers. Whether to actually recruit other machines, rather than just list them.'), ...
    arg({'verbose_output','VerboseOutput'},true, [], 'Verbose output. Whether to display information about the acquisition process.'));


if isempty(hostnames) || isequal(hostnames,{''})
    hostnames = par_globalpool; end

if isempty(hostnames)
    disp('The list of hostnames to connect to is empty; exiting.');
    return;
end

if ispc
    disp('Note: Acquiring workers from BCILAB will only work if your operating system provides an ''ssh'' command. From a Windows machine, you likely have to have some program installed to use this feature. Also note that your cluster needs to run Linux for this to work.'); 
    disp('You can always fall back to starting your worker processes (either as MATLAB instances or compiled binaries) by hand and listing their host:port''s in the worker pool (see Cluster Settings GUI).'); 
end

if isempty(mcr_root)
    mcr_root = ''; end
if isempty(identity_file)
    identity_file = ''; end

% sanity checks
if ~iscellstr(hostnames)
    error('The hostnames must be given as a cell array of strings.'); end
if ~isscalar(processors_per_node) || ~isnumeric(processors_per_node)
    error('The processors_per_node parameter must be given as a number.'); end
if ~islogical(share_resources)
    error('The share_resources parameter must be either true or false.'); end
if ~isscalar(sharing_port) || ~isnumeric(sharing_port)
    error('The sharing_port parameter must be given as a number.'); end
if ~ischar(logging_path)
    error('The logging_path parameter must be given as a string.'); end
if ~isempty(startup_prefix) && ~ischar(startup_prefix)
    error('The startup_prefix parameter must be given as a string.'); end
if ~ischar(startup_command)
    error('The startup_command parameter must be given as a string.'); end
if ~islogical(binary_worker)
    error('The binary_worker flag must be either true or false.'); end
if ~ischar(mcr_root)
    error('The mcr_root parameter must be given as a string.'); end
if ~ischar(identity_file)
    error('The identity_file parameter must be given as a string.'); end
if ~isscalar(max_cpuload) || ~isnumeric(max_cpuload)
    error('The max_cpuload parameter must be given as a number.'); end
if ~isscalar(min_memory) || ~isnumeric(min_memory)
    error('The min_memory parameter must be given as a number.'); end
if ~isscalar(load_window) || ~isnumeric(load_window)
    error('The load_window parameter must be given as a number.'); end
if ~isscalar(shutdown_timeout) || ~isnumeric(shutdown_timeout)
    error('The shutdown_timeout parameter must be given as a number.'); end
if ~isempty(identity_file)
    identity_file = ['-i ' identity_file]; end
window_remap = [1 1 1 2 2 2 2 2 2 2 3 3 3 3 3];
load_window = window_remap(min(load_window,15));

% generate in the bcilab-specific startup command
if strcmp(startup_command,'autogenerate')
    if binary_worker
        startup_command = sprintf('cd %s',env_translatepath('bcilab:/build/distrib')); 
    else
        startup_command = sprintf('%s; cd %s; bcilab %s worker {%s,1,%d} parallel {}',startup_prefix, env_translatepath('bcilab:/'), tracking.configscript,'%d',shutdown_timeout); 
    end
end

% use a port that depends on the startup command for sharing resources
if sharing_port == 0
    sharing_port = 10000 + mod(23457+hlp_fingerprint(startup_command),50000); end

% remove any port assignments
for i=1:length(hostnames)
    colons = hostnames{i}==':';
    if any(colons)
        hostnames{i} = hostnames{i}(1:find(colons)-1); end
end
% remove duplicates
hostnames = unique(hostnames);


% filter hostnames by machine availability, and retrieve machine stats
disp('Listing compute servers ...');
stats = [];
for host = hostnames(:)'    
    % collect system info
    [status,info] = system(sprintf('ssh %s -x %s "cat /proc/cpuinfo | grep ^processor | wc -l; cat /proc/loadavg; ps -A | grep ''%s'' | wc -l; cat /proc/meminfo | grep ^MemFree"',identity_file,host{1},avoid_proc));
    if ~status
        try
            % parse info
            lines = hlp_split(info,10);
            loadavg = hlp_split(lines{2},' ');
            meminfo = hlp_split(lines{4},': ');
            % store stats
            stats(end+1).hostname = host{1};
            try                
                stats(end).processors = str2num(lines{1});
            catch
                stats(end).processors = 4;
            end
            try
                stats(end).cpuload = str2num(loadavg{load_window});
            catch
                stats(end).cpuload = 0;
            end
            try
                stats(end).freemem = 1024*str2num(meminfo{2});
            catch
                stats(end).freemem = Inf;
            end
            try
                stats(end).avoidance = str2num(lines{3});
                if ~isscalar(stats(end).avoidance)
                    stats(end).avoidance = 0; end
            catch
                stats(end).avoidance = 0;
            end
        catch,end
    end
end

% start workers if necessary
pool = {};
if ~isempty(stats)
    fprintf('%.0f%% of requested hosts are available.\n',100*length(stats)/length(hostnames));
    
    % reorder nodes by increasing CPU load
    [dummy,idx] = sort([stats.cpuload]); %#ok<ASGLU>
    stats = stats(idx);
    hostnames = {stats.hostname};
    
    % start the workers...
    for k = 1:length(hostnames)
        host = hostnames{k};
        % check CPU load criterion
        if stats(k).cpuload > max_cpuload
            continue; end
        % check avoid criterion
        if stats(k).avoidance
            continue; end
        % determine number of ports (=procs) to consider
        if processors_per_node == 0
            % as many as there are cores
            rangelen = ceil(stats(k).processors/matlabthreads);
        else
            % at most 4x as many as (#cores/#matlabthreads)
            if ~isempty(stats(k).processors) && stats(k).processors > 0
                rangelen = min(ceil(stats(k).processors/matlabthreads)*4,processors_per_node);
            else
                rangelen = processors_per_node;
            end
        end
        % determine starting port range to consider
        if share_resources
            % the worker process may be shared between users
            rangestart = sharing_port;
        else
            while true
                % chose a random start port that does not overlap with the share port in a range of
                % 512 successive ports
                rangestart = 10000+50000*rand();
                if isempty(intersect(rangestart:rangestart+512-1,sharing_port:sharing_port+512-1))
                    break; end
            end
        end
        
        % try to start them...
        for port=rangestart:rangestart+rangelen-1
            if verbose_output
                fprintf('  acquiring worker at %s on port %d...\n', host, port); end
            % construct command
            if isempty(logging_path)
                logportion = '';
            else
                try
                    % ensure that the logging path exists
                    io_mkdirs([logging_path filesep],{'+w','a'});
                catch
                    disp_once(['Warning: the logging path "' logging_path '" does not exist locally and could not be created.\nPlease make sure that it exists on all worker machines, as otherwise the computation will not start.']);
                end
                logportion = sprintf(' > %s/%s_%d.out',logging_path,host,port);
            end
            if binary_worker
                if isempty(mcr_root)
                    % search it 
                    [major,minor] = mcrversion;
                    mcr_tail = sprintf([filesep 'MATLAB_Compiler_Runtime' filesep 'v%d%d'],major,minor);
                    possible_locations = {'/opt/MATLAB','/usr/local/MATLAB','/usr/common/MATLAB',hlp_homedir,[hlp_homedir filesep 'MATLAB'],'C:\\Program Files\\MATLAB','C:\\Program Files (x86)\\MATLAB','/Applications/MATLAB/'};
                    for p=1:length(possible_locations)
                        if exist([possible_locations{p} mcr_tail],'dir')
                            mcr_root = [possible_locations{p} mcr_tail]; end
                    end
                    if isempty(mcr_root)
                        fprintf('MATLAB compiler runtime v%d%d not found; if you get an error subsequently, please make sure that it is installed in a recognized location (or pass the mcr_root as an acquire options).\n',major,minor); end
                end
                if isempty(mcr_root)
                    % assume that the binary finds all necessary MCR files
                    command = ['ssh ' identity_file ' -x ' host '  "' startup_command '; ./' binary_name ' worker ''{' num2str(port) ',1,' num2str(shutdown_timeout) '}'' parallel {} ' logportion '" &'];
                else
                    % use the run script to point to the correct MCR install location
                    command = ['ssh ' identity_file ' -x ' host ' "' startup_command '; ./run_' binary_name '.sh ' mcr_root ' worker ''{' num2str(port) ',1,' num2str(shutdown_timeout) '}'' parallel {} ' logportion '" &'];
                end
            else
                % run the MATLAB startup command (and substitute the port into it)
                command = sprintf('ssh %s -x %s "%s -nodisplay -r ''%s''%s" &', ...
                    identity_file, host, matlab_command, sprintf(startup_command,port),logportion);
            end
            % issue startup command, if necessary
            if recruit_machines && port_free(host,port)
                % port free: check for memory criterion
                if stats(k).freemem > min_memory
                    % start, and add to pool
                    system(command);
                    stats(k).freemem = stats(k).freemem - min_memory;
                    pool{end+1} = sprintf('%s:%d',host,port);
                end
            else
                % worker is already running: add to pool
                pool{end+1} = sprintf('%s:%d',host,port);
            end
        end
    end    
    fprintf('%i processor slots acquired.\n',length(pool));
else
    disp('None of the listed machines is available; please make sure that:');
    disp(' * you can ssh into these machines without a password prompt"'); 
    disp(' * the files /proc/cpuinfo, /proc/meminfo, and /proc/loadavg are available on them"'); 
    disp(' * the commands ps, cat, grep, and wc are working on them"'); 
end

if nargout == 0
    % assign to global variable, if necessary...
    tracking.parallel.pool = pool(:)';
end


% check whether the given port on the specified host is free
function result = port_free(host,port)
try
    java.net.Socket(host,port);
    result = false;
catch
    result = true;
end
