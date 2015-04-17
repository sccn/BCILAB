function [harvested_addresses,logpaths] = par_getworkers_qsub(varargin)
% Acquire workers on some remote machines and return hostnames and ports of those that are available.
% Pool = par_getworkers_qsub(...)
%
% This function attempts to start the desired number of workers on a Linux cluster using the qsub
% command (also supports some other job managers). See aksi 
%
% In:
%   --- resources to acquire ---
%   NumWorkers : Number of worker processes to launch. (default: 8)
%
%   MatlabThreads : Number of MATLAB threads. This is the number of threads that each worker uses internally. (default: 4)
%
%   --- qsub-specific options ---
%   SubmitNode : Machine which can submit jobs. If nonempty, this function will attempt to ssh into 
%                the given node in order to submit a job. (default: '')
%
%   JobManager : Job manager system. This is the type of job manager to submit to. (default: 'SGE')
%
%   Queues : Queues to submit to. This is a cell array of strings that lists the queues to which jobs
%            shall be submitted to (in a round-robin manner or in a merged manner). If empty, no
%            queue is specified. (default: {})
%
%   HarvestTimeout : Harvesting timeout. Any workers that have not fully launched within this
%                    timeout will be abandoned (and will terminate themselves). Note that when
%                    multiple MATLABs try to launch simultaneously, it can take longer than normal.
%                    (default: 300)
%
%   --- startup commands for MATLAB ---
%   MATLABCommand : Command to start MATLAB. This is the command to run on the command line to start 
%                   the desired version of MATLAB (does not include the script launch). If this is
%                   set to ''autogenerate'' the same path that runs the par_getworkers_qsub command
%                   will be used (useful on clusters with identical file systems). A good fallback
%                   on most installations is to set it to ''matlab''. (default: 'autogenerate')
%
%   StartupPrefix : Startup prefix. Any MATLAB startup lines to run before running the main BCILAB
%                   command. (default: '')
%
%   StartupCommand : BCILAB startup command. If the binary_worker parameter (Use compiled workers) 
%                    is true, this is a sequence of shell commands, otherwise it is a sequence of
%                    MATLAB statements which contain the placeholder %d where the TCP port # would
%                    be inserted; If set to ''autogenerate'', a startup command specific to the
%                    current BCILAB environment will be generated. (default: 'autogenerate')
%
%   NoDisplay : Start MATLAB without display. (default: true)
%
%   CleanMATLABPath : Run worker from a clean MATLAB path. (default: true)
%
%   --- binary commands for compiled workers ----
%	BinaryWorker : Use compiled workers. Whether to start a given binary worker implementation
%                  rather than the MATLAB worker. (default: false)
%
%   MatlabCompilerRoot : MATLAB Compiler directory. Installation path of the MATLAB compiler runtime; 
%                        if set to autogenerate, a few locations will be searched for an MCR that
%                        corresponds to the MATLAB version that is running this function (this
%                        assumes that the compiler is installed in the same directory remotely as it
%                        is locally). If nothing is found, it is assumed that the MCR has been
%                        installed into the system path and is found automatically. (default:
%                        autogenerate)
%
%   BinaryName : Name of the binary. The name of the BCILAB binary (if using compiled workers).
%                (default: 'build')
%
%   --- arguments for the workers ---
% 	StartPort : Start of the worker port range. If the port is already in use, the next free one
%               will be chosen, until the permitted portrange is exceeded. if specified as 0, a free
%               port is chosen directly. (default: 23547)
%
%   NumPorts : Number of available ports. Number of successive ports to try, including the start
%              port; if 1, only the desired port will be used. If 0, the port range will be set to
%              the number of cores in the machine divided by the worker''s matlab threads parameter.
%              (default: 0)
%
%   ShutdownTimeout : Shutdown timeout. If non-zero, the worker will be shut down if it has not
%                     received a heartbeat signal from a client in the given time frame, in seconds.
%                     For this to work, the function env_acquire_cluster should be used as it sets
%                     up the heartbeat timer. (default: 300)
%
%   ---- worker tracking and logging ---
%   JobIdFormat : Job ID format. This is a pattern according to which job ids are generated. 
%                 (default: 'jobmanager_%user_%host_b%batch_i%num')
%
%   LoggingPath : Logging path pattern. For each worker job a new file will be created according to
%                 this pattern. (default: 'home:/.bcilab/logs/workers/qsub-%jobid.log')
%
%
% Out:
%   Pool : cell array of 'hostname:port' strings specifying the list of available machines
%          Note: in the case that no return value is requested, the global variable 
%                tracking.parallel.pool will receive this result. This is the recommended way to
%                use par_getworkers_qsub, as par_schedule uses this pool by default.
%
%  Logpaths : cell array of file paths of the logfiles corresponding to the workers in pool
%
% See also:
%   par_worker
%
% Examples:
%   % launch three workers on any of the queues q1 to q8, using computing as the submit host
%   par_getworkers_qsub('NumWorkers',3,'Queues',{'q1','q2','q3','q4','q5','q6','q7','q8'},'SubmitNode','computing')
%
% Notes:
%   Some details of this function were inspired by the qsub module of FieldTrip.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-06-27

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
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
opts = arg_define(varargin,...
    ... % resources to acquire
    arg({'num_workers','NumWorkers'},8,[1 1 32 1024],'Number of worker processes to launch.'), ...
    arg({'matlab_threads','MatlabThreads'},4,uint32([1 1 8 32]),'Number of MATLAB threads. This is the number of threads that each worker uses internally.'), ...
    ... % qsub-specific options
    arg({'submit_node','SubmitNode'},'',[],'Machine which can submit jobs. If nonempty, this function will attempt to ssh into the given node in order to submit a job.'), ...
    arg({'job_manager','JobManager'},'SGE',{'SGE','OGE','Torque','PBS','SLURM','LSF','OwnMachine'},'Job manager system. This is the type of job manager to submit to.'), ...
    arg({'queues','Queues'},{},[],'Queues to submit to. This is a cell array of strings that lists the queues to which jobs shall be submitted to (in a round-robin manner or in a merged manner). If empty, no queue is specified.','type','cellstr','shape','row'), ...
    arg({'harvest_timeout','HarvestTimeout'},300,[],'Harvesting timeout. Any workers that have not fully launched within this timeout will be abandoned (and will terminate themselves). Note that when multiple MATLABs try to launch simultaneously, it can take longer than normal.'), ...
    arg({'harvest_ips','HarvestIPs'},true,[],'Harvest IP addresses instead of hostnames. Whether to harvest and return the IP addresses of the workers rather than their hostnames.'), ...
    ... % startup commands for MATLAB
    arg({'matlab_command','MATLABCommand'},'autogenerate', [], 'Command to start MATLAB. This is the command to run on the command line to start the desired version of MATLAB (does not include the script launch). If this is set to ''autogenerate'' the same path that runs the par_getworkers_qsub command will be used (useful on clusters with identical file systems). A good fallback on most installations is to set it to ''matlab''.'), ...
    arg({'startup_prefix','StartupPrefix'},'',[], 'Startup prefix. Any MATLAB startup lines to run before running the main BCILAB command.'),...
    arg({'startup_command','StartupCommand'},'autogenerate',[], 'BCILAB startup command. If the binary_worker parameter (Use compiled workers) is true, this is a sequence of shell commands, otherwise it is a sequence of MATLAB statements which contain the placeholder %d where the TCP port # would be inserted; If set to ''autogenerate'', a startup command specific to the current BCILAB environment will be generated.'),...
    arg({'no_display','NoDisplay'},true,[], 'Start MATLAB without display.','guru',true),...
    arg({'clean_path','CleanMATLABPath'},true, [], 'Run worker from a clean MATLAB path.'), ...
    ... % binary commands for compiled workers
    arg({'binary_worker','BinaryWorker'},false,[], 'Use compiled workers. Whether to start a given binary worker implementation rather than the MATLAB worker.'), ...
    arg({'mcr_root','MatlabCompilerRoot'},'autogenerate', [], 'MATLAB Compiler directory. Installation path of the MATLAB compiler runtime; if set to autogenerate, a few locations will be searched for an MCR that corresponds to the MATLAB version that is running this function (this assumes that the compiler is installed in the same directory remotely as it is locally). If nothing is found, it is assumed that the MCR has been installed into the system path and is found automatically.'), ...
    arg({'binary_name','BinaryName'},'build', [], 'Name of the binary. The name of the BCILAB binary (if using compiled workers).'), ...
    ... % arguments for the workers
    arg({'start_port','StartPort'}, 23547, [0 1025 65535 65535], 'Start of the worker port range. If the port is already in use, the next free one will be chosen, until the permitted portrange is exceeded. if specified as 0, a free port is chosen directly.'), ...
    arg({'num_ports','NumPorts'}, 0, [0 0 64 65536], 'Number of available ports. Number of successive ports to try, including the start port; if 1, only the desired port will be used. If 0, the port range will be set to the number of cores in the machine divided by the worker''s matlab threads parameter.'), ...
    arg({'shutdown_timeout','ShutdownTimeout'},300,[0 60 3600 Inf],'Shutdown timeout. If non-zero, the worker will be shut down if it has not received a heartbeat signal from a client in the given time frame, in seconds. For this to work, the function env_acquire_cluster should be used as it sets up the heartbeat timer.'), ...
    ... % worker tracking and logging
    arg({'jobid_format','JobIdFormat'},'jobmanager_%user_%host_b%batch_i%num', [], 'Job ID format. This is a pattern according to which job ids are generated.','guru',true), ...
    arg({'logging_path','LoggingPath'},'home:/.bcilab/logs/workers/qsub_b%batch/%jobid.log', [], 'Logging path pattern. For each worker job a new file will be created according to this pattern.'));

% no arguments were passed? bring up GUI dialog
if isempty(varargin)
    opts = arg_guidialog;    
    if isempty(opts)
        return; end % -> user clicked cancel
end

% copy options to workspace
[num_workers,matlab_threads,submit_node,job_manager,queues,harvest_timeout,harvest_ips,matlab_command,startup_prefix,startup_command,no_display,clean_path,binary_worker,mcr_root,binary_name,start_port,num_ports,shutdown_timeout,jobid_format,logging_path] = arg_toworkspace(opts);

% pre-generate job ids according to the jobid_format
if hlp_matlab_version >= 712
    rng shuffle
else
    rand('seed',sum(100*clock)); %#ok<RAND>
end
batchid = num2str(10000 + round(rand*89999));
username = char(java.lang.System.getProperty('user.name'));
hostname = hlp_hostname;
job_ids = cell(1,num_workers);
for k=1:num_workers
    job_ids{k} = sanitize_name(strrep_multi(jobid_format,'%user',username,'%host',hostname,'%batch',batchid,'%num',num2str(k,'%03i'))); end


% process the matlab command
if isempty(matlab_command) && ~binary_worker
    error('The given MATLAB command (%s) must not be empty. Note that you can set it to ''autogenerate'' if unsure.'); end
if strcmp(matlab_command,'autogenerate') 
    % use the currently running MATLAB executable
    matlab_command = [matlabroot filesep 'bin' filesep 'matlab']; 
    if ~exist(matlab_command,'file')
        fprintf('WARNING: the MATLAB binary (%s) does not seem to exist on the local file system. Note that you can override it in the options.\n',matlab_comand); end
else
    % sanity-check the user-provided MATLAB command (trick by R. Oostenveld)
    if ~ispc && system(sprintf('which %s > /dev/null', matlabcmd))==1
        fprintf('WARNING: the given MATLAB command (%s) does not seem to exist on the local system; make sure that it exists on the remote cluster.\n',matlab_command); end
end


% process the startup prefix
if clean_path && ~binary_worker
    startup_prefix = ['restoredefaultpath; ' startup_prefix]; end


% process the bcilab-specific startup command
if strcmp(startup_command,'autogenerate')
    if binary_worker
        startup_command = sprintf('cd %s',env_translatepath('bcilab:/build/distrib'));
    else
        startup_command = sprintf('%s; cd %s; bcilab %s worker {%d,%d,%d} parallel {}',startup_prefix, env_translatepath('bcilab:/'), tracking.configscript, start_port, num_ports, shutdown_timeout);
    end
end


% process the MCR root path
if binary_worker
    if strcmp(mcr_root,'autogenerate')
        [major,minor] = mcrversion;
        mcr_tail = sprintf([filesep 'MATLAB_Compiler_Runtime' filesep 'v%d%d'],major,minor);
        possible_locations = {'/opt/MATLAB','/usr/local/MATLAB','/usr/common/MATLAB',hlp_homedir,[hlp_homedir filesep 'MATLAB'],'C:\\Program Files\\MATLAB','C:\\Program Files (x86)\\MATLAB','/Applications/MATLAB/'};
        for p=1:length(possible_locations)
            if exist([possible_locations{p} mcr_tail],'dir')
                mcr_root = [possible_locations{p} mcr_tail]; end
        end
        if isempty(mcr_root)
            fprintf('MATLAB compiler runtime v%d%d not found; if you get an error subsequently, please make sure that it is installed in a recognized location (or pass the mcr_root as an acquire options).\n',major,minor); end
    elseif ~isempty(mcr_root) && ~exist(mcr_root,'dir')
        fprintf('NOTE: The matlab compiler root (%s) does not exist on the local file system; make sure that it is accessible to remote workers.\n',mcr_root);
    end
end


% process the binary name
if binary_worker && isempty(mcr_root)
    binary_name = ['run_' binary_name '.sh']; end


% process the logging paths
if isempty(logging_path)
    error('The logging path must not be empty (since log files will be used to identify what workers launched successfully).'); end
logpaths = cell(1,length(job_ids));
for k=1:length(job_ids)
    logpaths{k} = env_translatepath(strrep_multi(logging_path,'%jobid',job_ids{k},'%batch',batchid,'%user',username,'%host',hostname,'%num',num2str(k,'%03i')));
    try
        % ensure that the logging path exists
        io_mkdirs(logpaths{k},{'+w','a'});
    catch
        fprintf('NOTE: the logging path "%s" could not be created.\nPlease make sure that it exists on all worker machines, as otherwise the computation will not start.\n',tmppath);
    end    
end


% build MATLAB command-line options
matlab_options = ' -logfile %logpath ';
if no_display
    matlab_options = [matlab_options ' -nodisplay ']; end
if matlab_threads == 1 && hlp_matlab_version > 708
    matlab_options = [matlab_options ' -singleCompThread ']; end


% build the job command line (different for binary vs. interpreted)
if binary_worker
    job_cmdline = [startup_command '; ./' binary_name mcr_root ' worker ''{' num2str(start_port) ',' num2str(num_ports) ',' num2str(shutdown_timeout) '}'' parallel {} > ' logging_path];
else
    job_cmdline = [matlab_command ' ' matlab_options ' -r ''' startup_command ''''];
end

% build qsub options
qsub_options = '';

% build the launch command (note: only SGE has been tested so far)
switch job_manager
    case {'SGE','OGE'}
        if ~isempty(queues)
            qsub_options = [qsub_options ' -q %queue ']; end
        launch_command = ['echo "' job_cmdline '" | qsub -N %jobid ' qsub_options ' -cwd -o %outputdir -e %outputdir'];
    case {'Torque','PBS'}
        if ~isempty(queues)
            qsub_options = [qsub_options ' -q %queue ']; end
        launch_command = ['echo "' job_cmdline '" | qsub -N %jobid ' qsub_options ' -d "%outputdir" -o "%outputdir" -e "%outputdir"'];
    case 'LSF'
        if ~isempty(queues)
            qsub_options = [qsub_options ' -q %queue ']; end
        launch_command = ['echo "' job_cmdline '" | bsub -J %jobid ' qsub_options ' -o %outputdir%jobid.out -e %outputdir%jobid.err'];
    case 'SLURM'
        if ~isempty(queues)
            qsub_options = [qsub_options ' --partition=%queue ']; end
        launch_command = ['srun --job-name=%jobid ' qsub_options ' --output=%outputdir%jobid.out --error=%outputdir%jobid.err ' job_cmdline];
    case 'OwnMachine'
        launch_command = job_cmdline;
    otherwise
        error('Unrecognized job manager type: %s',job_manager);
end


% merge multiple queues if applicable
if length(queues) > 1 && any(strcmp(job_manager,{'SGE','OGE'}))
    queues = sprintf('%s,',queues{:});
    queues = {queues(1:end-1)};
end

io_mkdirs(env_translatepath('temp:/job_management/'),{'+w','a'});
output_dir = env_translatepath('temp:/job_management');
if any(output_dir== ' ') 
    fprintf('Note: your BCILAB temp directory contains a space character; job output will be dumped into your home directory instead.\n');
    output_dir = hlp_homedir;
end
output_dir = [output_dir filesep];
queue = '';

% submit the jobs
for k=1:num_workers
    % determine per-job parameters
    job_id = job_ids{k};
    if ~isempty(queues)
        queue = queues{1+mod(k,length(queues))}; end
    
    % perform substitutions in the launch command
    issue_commandline = strrep_multi(launch_command, ...
        '%jobid',job_id, '%outputdir',output_dir, '%queue',queue, '%logpath',logpaths{k});
    
    % invoke command and check for errors
    fprintf('Scheduling worker #%i (%s): %s...\n',k,job_id,issue_commandline);
    if isempty(submit_node)
        % regular case: invoke it locally
        [status, output] = system(issue_commandline);
        if status
            error(output); end
    else
        % we need to invoke the command on a remote machine to submit the job
        submit_file = env_translatepath(['temp:/job_management/' job_id '.sh']);
        fid = fopen(submit_file,'w');
        if fid==-1
            error('Could not open submit file %s for writing; please check your permissions.',submit_file); end
        fprintf(fid,issue_commandline);
        fclose(fid);
        system(['chmod +x ' submit_file]);
        [status, output] = system(['ssh -x ' submit_node ' ' submit_file]);
        if status
            error(output); end
    end        
end

fprintf('\nWaiting for workers to start up to establish connections...\n');
harvested_addresses = par_parse_logfiles(logpaths, harvest_timeout, harvest_ips);
fprintf('Launched %i workers: %s\n',length(harvested_addresses),hlp_tostring(harvested_addresses));
    
if nargout == 0
    % assign to global variable, if necessary...
    par_globalsetting('pool',harvested_addresses(:)');
    par_globalsetting('logfiles',logpaths(:)');
end

function str = strrep_multi(str,varargin)
for k=1:2:length(varargin)
    str = strrep(str,varargin{k},varargin{k+1}); end

function s = sanitize_name(s)
s(~((s>='a'&s<='z') | (s>='A'&s<='Z') | (s>='0'&s<='9') | s=='-')) = '_';
