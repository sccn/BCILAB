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
opts = hlp_varargin2struct(varargin, ...
     'engine','global', ...
     'pool', 'global', ...
     'policy', 'global', ...
     'verbosity', 'global', ...
     'receiver_backlog', 5, ...
     'receiver_timeout', 1, ...
     'reschedule_interval',5, ...
     'keep', false, ...
     'pushscope', true ...
    );

% look up global settings, if requested
if strcmp(opts.engine,'global')
    opts.engine = par_globalsetting('engine'); end
if strcmp(opts.pool,'global')
    opts.pool = par_globalsetting('pool'); end
if strcmp(opts.policy,'global')
    opts.policy = par_globalsetting('policy'); end
if strcmp(opts.verbosity,'global')
    opts.verbosity = par_globalsetting('verbosity'); end
if isa(opts.policy,'function_handle')
    opts.policy = char(opts.policy); end
if strcmp(opts.engine,'BLS') && isempty(opts.pool)
    opts.engine = 'local'; end

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
        error('Unsupported task format; please see documentation.');
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
    if opts.keep
        tmp = hlp_microcache('schedulers',@(varargin)Scheduler(varargin{:}),opts.pool,opts.policy,opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval),opts.verbosity,length(tasks));
    else
        tmp = Scheduler(opts.pool,opts.policy,opts.receiver_backlog,round(1000*opts.receiver_timeout),round(1000*opts.reschedule_interval),opts.verbosity,length(tasks));
    end
    sched = struct('sched',{tmp},'finisher',{onCleanup(@()tmp.clear())});
end

% serialize tasks for network transmission (and prepend order id)
if any(strcmp(opts.engine,{'BLS','Reference'}))
    for t=1:length(tasks)        
        tasks{t} = fast_encode(hlp_serialize([{t} tasks{t}])); end
end

% submit tasks for execution
switch opts.engine
    case 'local'
        sched = {};
        for t=1:length(tasks)
            sched(t) = {{t,tasks{t}{1}(tasks{t}{2:end})}}; end
    case 'ParallelComputingToolbox'
        sched = {};
        parfor t=1:length(tasks)
            sched(t) = {{t,tasks{t}{1}(tasks{t}{2:end})}}; end
    case 'BLS'
        % over the scheduler
        sched.sched.submit(tasks);
    case 'Reference'
        % evaluate locally, but go through the same evaluation function as the BLS workers
        for t=1:length(tasks)
            tasks{t} = fast_encode(par_evaluate(fast_decode(tasks{t}))); end
        % return the collected result in sched
        sched = struct('ReferenceResults',{tasks});
    otherwise
        error('Unsupported parallelization engine selected.');
end
