function [pool,logpaths] = par_getworkers(varargin)
% Acquire workers on some remote machines and return hostnames and ports of those that are available.
% Pool = par_getworkers(...)
%
% This function attempts to start the desired number of workers; offers multiple mechanisms to do
% so.
%
% In:
%   System : cell array of {Mechanism, Arguments...} where Mechanism is one of the following:
%            * 'ssh' : use ssh to launch workers on a list of given Linux machines
%            * 'qsub' : use qsub to submit jobs to a job manager such as Sun Grid Engine
%                       (also supports some non-qsub job managers)
%            the Arguments are name-value pairs accepted by the function par_getworkers_<Mechanism>
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
%   % use qsub to launch three workers on any of the queues q1 to q8, using computing as the submit host
%   par_getworkers({'qsub', 'NumWorkers',3, 'Queues',{'q1','q2','q3','q4','q5','q6','q7','q8'}, 'SubmitNode','computing')
%
%   % use ssh to launch 12 workers total on the three hosts computing-0-1, 0-2 and 0-3
%   par_getworkers({'ssh', 'Hostnames',{'computing-0-1','computing-0-2',computing-0-3'}, 'ProcessorsPerNode',4)
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

opts = arg_define(varargin,...
    arg_subswitch({'system','System'},{'qsub'},{'ssh',@par_getworkers_ssh,'qsub',@par_getworkers_qsub},'Job acquisition system to use. Different systems are available, including ssh (logging into each node and launching a number of workers) and qsub (submitting worker jobs to a job manager).'));

% dispatch to respective sub-function
[pool,logpaths] = feval(['par_getworkers_' opts.system.arg_selection],opts);
