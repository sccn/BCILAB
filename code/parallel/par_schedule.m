function [results,errors] = par_schedule(tasks,varargin)
% Schedule the given tasks across a pool of (possibly remote) workers.
% Results = par_schedule(Tasks, Options...)
%
% In:
%   Tasks : cell array of tasks; formatted as
%           * (evaluatable) string
%           * {function_handle, arg1,arg2,arg3, ...}
%           * struct('head',function_handle, 'parts',{{arg1,arg2,arg3, ...}})
%           (see also par_beginschedule for further details)
%
%   Options...: optional name-value pairs, same as in par_beginschedule, with the addition of
%
%               'scope': optional parallel scope; if this is a cell array of name-value pairs, 
%                        cluster resources will be acquired with these options for the duration of
%                        bci_train (and released thereafter). Options as in
%                        env_acquire_cluster.
%
% Out:
%   Results : cell array of results of the scheduled computations (evaluated tasks)
%   Errors  : cell array of exception structs for those results that could not be evaluated (in no particular order)
%
% See also:
%   par_worker, par_beginschedule, par_endschedule
%
% Notes:
%   Only the first output value of each task is taken and returned, though you can schedule 
%   hlp_wrapresults or hlp_getresult to get all or a specific output value of your task function.
%
% Example:
%   % run two computations (here as strings) in parallel on a pool of two nodes (assuming that MATLAB
%   % is running on those, executing the par_worker function)
%   results = par_schedule({'sin(randn(10))','exp(randn(10))'},'pool',{'localhost:32547','localhost:32548'})
% 
%   % as before, but pass the jobs as {function,arguments...}
%   results = par_schedule({@sin,randn(10)},{@exp,randn(10)}},'pool',{'localhost:32547','localhost:32548'})
%   
%   % as before, but do not destroy and re-create the scheduler between calls to par_schedule
%   results = par_schedule({@sin,randn(10)},{@exp,randn(10)}},'pool',{'localhost:32547','localhost:32548'},'keep'true)
%
%                       Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                       2010-08-29

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

% optionally set up a scoped parallel run on a cluster
opts = hlp_varargin2struct(varargin,'scope',[]);
if iscell(opts.scope)
    if env_acquire_cluster(opts.scope{:})
        releaser = onCleanup(@()env_release_cluster); end
end

% if there's only one task we run locally...
if length(tasks) == 1
    varargin = [varargin {'engine','local'}]; end
    
% schedule tasks...
id = par_beginschedule(tasks,varargin{:});
[results,errors] = par_endschedule(id,varargin{:});

if ~isempty(errors) && nargout <= 1
    rethrow(errors{1}{2}); end
