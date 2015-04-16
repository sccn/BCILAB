function sched = par_beginschedule(tasks,varargin)
% Begin the scheduling of some set of tasks across a pool of (possibly remote) workers.
% Id = par_beginschedule(Tasks,Options...)
%
% Returns a scheduler handle to wait for and obtain results upon completion.
%
% In:
%   Tasks : cell array of tasks; each cell can be formatted as 
%           * (evaluatable) string
%           * {function_handle, arg1,arg2,arg3, ...}
%           * struct('head',function_handle, 'parts',{{arg1,arg2,arg3, ...}})
%
%   Options...: optional name-value pairs, with possible names:
%               'engine': parallelization engine to use, can be one of:
%                         'global': select the global setting (tracking.parallel.engine)
%                         'local': do all computations locally, skipping serialization
%                         'BLS': use the BCILAB Scheduler (uses the resources specified in the pool argument) (default)
%                         'Reference': local reference implementation for testing BLS (using the same task serialization mechanism)
%                         'ParallelComputingToolbox': use the Mathworks Parallel Computing Toolbox (tm); uses resources allocated via the matlabpool command or a configuration file
%
%               'pool': pool of workers to consider for the BLS scheduler (default: 'global')
%                       if 'global', the global setting (tracking.parallel.pool) will be chosen
%                       (with the BLS engine, an empty pool implies local computation)
%
%               'policy': name of the scheduling policy function for the BLS (default: 'global')
%                         if 'global', the global setting (tracking.parallel.policy) will be chosen
%
%               'pushscope' : whether to "push" the current symbol scope (see hlp_scope and hlp_resolve)
%                             over the network (default: true)
%
%               'keep': keep this scheduler alive for later re-use (default: false)
%                       if false, the scheduler will be destroyed after use, and re-created during the next run
%
% Out:
%   Id : output handle; used to collect the results
%
% See also:
%   par_endschedule, par_worker, par_schedule
%
% Example:
%   % schedule two computations across a pool of some IP:port addresses (assuming that MATLAB is running there,
%   % and is executing the par_worker() function
%   id = par_beginschedule({'sin(randn(10))','exp(randn(10))'}, 'pool',{'192.168.1.1:23547','192.168.1.2:23547','192.168.1.2:23548','192.168.1.3:23547'});
%   ... optionally do something in the meantime
%   results = par_endschedule(id);
%
% Expert note:
%  The 'keep' option is easy to use with the wrapper function par_schedule; otherwise, the following holds:
%   * if passed as true to par_beginschedule, it *must* also be passed as true to par_endschedule
%   * nested schedules are not allowed if they use the same worker pool
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-08-29

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

% read options
opts = arg_define('allow-unlisted-names',varargin, ...
    arg_norep({'tasks','Tasks'},[],[],'Tasks to execute.','type','expression'), ...
    arg({'engine','ParallelEngine'},'global',{'global','local','BLS','Reference','ParallelComputingToolbox'}, 'Parallel engine to use. This can either be one of the supported parallel engines (BLS for BCILAB Scheduler, Reference for a local reference implementation, and ParallelComputingToolbox for a PCT-based implementation), or local to skip parallelization altogether, or global to select the currently globally selected setting (in the global tracking variable).'), ...
    arg({'pool','WorkerPool'},'global',[], 'Worker pool to use. This is typically a cell array, but can also be the string ''gobal'', which stands for the currently globally set up worker pool (see global tracking variable).','type','expression'), ...
    arg({'policy','ReschedulingPolicy'},'global',[], 'Rescheduling policy. This is the name of the rescheduling policy function that controls if and when tasks are being rescheduled. If set to global, the current global setting will be used.'), ...
    arg({'verbosity','Verbosity'},'global',{'global','0','1','2','3'}, 'Verbosity level. Verbosity level for the scheduling process, or global to use the global setting.'), ...
    arg({'receiver_backlog','ReceiverBacklog'},5,uint32([0 1 7 10]), 'Max backlogged connections at receiver. Maximum number of connections that can be in the backlog of the result receiver. '), ...
    arg({'receiver_timeout','ReceiverTimeout'},1,[0 1 60 3600], 'Receiver timeout. In seconds. When a worker that has previously announced that it is ready to transmit results does not start transmitting results within this time period after pinged, the connection will be discarded and the worker has to re-connect.'), ...    
    arg({'reschedule_interval','RescheduleInterval'},5,[0 1 15 3600], 'Rescheduling interval. In seconds. Time between periodic checks if jobs have to be re-scheduled.'), ...
    arg({'keep','RetainScheduler'},false,[], 'Retain scheduler instance. If true, the scheduler instance will be reused on the next call to par_beginschedule. Can be more efficient, but does not allow for nested calls to par_beginschedule.'), ...
    arg({'pushscope','PushScope'},true,[], 'Push variable scope. If true, the current variable scope (established via hlp_scope) will be pushed over the network; necessary for some features to run seamlessly.'));

% look up global settings, if requested
if strcmp(opts.engine,'global')
    opts.engine = par_globalsetting('engine'); end
if strcmp(opts.pool,'global')
    opts.pool = par_globalsetting('pool'); end
if strcmp(opts.policy,'global')
    opts.policy = par_globalsetting('policy'); end
if strcmp(opts.verbosity,'global')
    opts.verbosity = par_globalsetting('verbosity'); end
if ischar(opts.verbosity)
    opts.verbosity = str2num(opts.verbosity); end
if isa(opts.policy,'function_handle')
    opts.policy = char(opts.policy); end
if strcmp(opts.engine,'BLS') && isempty(opts.pool)
    opts.engine = 'local'; end

% validate options
if ~iscellstr(opts.pool)
    error('The given worker pool must be a cell array of strings, but was: %s',hlp_tostring(opts.pool)); end
if ~ischar(opts.policy)
    error('The given rescheduling policy must be a string, but was: %s',hlp_tostring(opts.policy)); end
if ~exist(opts.policy,'file')
    error('The given rescheduling policy refers to a non-existing file: %s',opts.policy); end

% canonlicalize task format to {function-handle,arguments...}
for t=1:length(tasks)
    if isfield(tasks{t},{'head','parts'})
        % task given as Mathematica-style expression struct (see expressions/exp_*)
        tasks{t} = [{tasks{t}.head} tasks{t}.parts];
    elseif ischar(tasks{t})
        % task given as string
        tasks{t} = {@eval,tasks{t}};
    elseif ~(iscell(tasks{t}) && ~isempty(tasks{t}) && isa(tasks{t}{1},'function_handle'))
        % incorrect task format...
        error('Unsupported task format, please see documentation: %s',hlp_tostring(tasks{t},10000));
    end
end

% push current symbol context over the network
if opts.pushscope && ~strcmp(opts.engine,'local')
    % get the current scope
    scope = hlp_resolveall;
    % and wrap a hlp_scope() around the task
    for t=1:length(tasks)
        tasks{t} = [{@hlp_scope, scope} tasks{t}]; end
end

% create a scheduler (Java code, see dependencies/Scheduling-*)
if strcmp(opts.engine,'BLS')
    if opts.verbosity>0 %#ok<*ST2NM>
        fprintf('Creating scheduler...\n'); end
    if opts.keep
        tmp = hlp_microcache('schedulers',@(varargin)Scheduler(varargin{:}),opts.pool,opts.policy,'par_accept_results',opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval),opts.verbosity,length(tasks));
    else
        tmp = Scheduler(opts.pool,opts.policy,'par_accept_results',opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval),opts.verbosity,length(tasks));
    end
    sched = struct('sched',{tmp},'finisher',{onCleanup(@()tmp.clear())});
end

% serialize tasks for network transmission (and prepend order id)
if any(strcmp(opts.engine,{'BLS','Reference'}))
    if opts.verbosity>0
        t0=tic; fprintf('Serializing %d tasks...\n', length(tasks)); end
    for t=1:length(tasks)        
        tasks{t} = fast_encode(hlp_serialize([{t} tasks{t}])); end
    if opts.verbosity>0
        fprintf('Tasks serialized (%.1f seconds).\n', toc(t0)); end
end

% submit tasks for execution
switch opts.engine
    case 'local'
        sched = cell(1,length(tasks));
        for t=1:length(tasks)
            sched(t) = {{t,tasks{t}{1}(tasks{t}{2:end})}}; end
    case 'ParallelComputingToolbox'
        sched = {};
        parfor t=1:length(tasks)
            sched(t) = {{t,tasks{t}{1}(tasks{t}{2:end})}}; end
    case 'BLS'
        % over the scheduler
        if opts.verbosity>0
            fprintf('Submitting tasks to scheduler...\n'); end
        sched.sched.submit(tasks);
    case 'Reference'
        % evaluate locally, but go through the same evaluation function as the BLS workers
        for t=1:length(tasks)
            tasks{t} = fast_encode(par_evaluate(fast_decode(tasks{t}))); end
        % return the collected result in sched
        sched = struct('ReferenceResults',{tasks});
    otherwise
        error('Unsupported parallelization engine selected: %s',hlp_tostring(opts.engine));
end
