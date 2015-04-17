function res = par_globalsetting(name,val)
% Set or get a global setting for parallel task scheduling.
% Result = par_globalsetting(SettingName,Value)
%
% In:
%   SettingName : name of the setting to look up or write; the most relevant settings are:
%                 'engine' : the compute engine to use; can be one of:
%                             'local': run within the local MATLAB context, in a serial manner 
%                                      (using a for loop) (default)
%                             'BLS': dispatch jobs to a collection of known workers, using the 
%                                    BCILAB Scheduler (BLS)
%                             'Reference': run locally in a for loop, but use the same task 
%                                          serialization and evaluation mechanism used by BLS, for
%                                          debugging
%                             'ParallelComputingToolbox': run using the parallel computing toolbox
%                                                         (using a parfor loop); this requires that 
%                                                         the toolbox is correctly configured
%
%                 'pool'   : pool of known worker machines for use with the BLS scheduler; 
%                            cell array of {'host:port','host:port',...} (default: {})
%
%                 'policy' : the current reschedule policy, for use with the BLS scheduler;
%                            (default: par_reschedule_policy)
%
%                 'verbosity' : the current verbosity level of the parallel computing engine
%                               0=only errors, 1=verbose, 2=extremely verbose (default: 0)
%
%                 'logfiles' : logfiles of the known worker processes (default: {})
%
%   Value : new value if the given setting shall be overridden (otherwise omitted)
%
% Out:
%   Result : previous value of the given setting
%
% Notes:
%   This function reads or writes contents of the global variable tracking.parallel. The preferred 
%   way to access these global variables is this function.
%
% See also:
%   par_schedule, par_worker
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2013-10-30

% Copyright (C) Christian Kothe, SCCN, 2013, christian@sccn.ucsd.edu
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


global tracking

if ~any(strcmp(name,{'engine','pool','policy','verbosity','logfiles'}))
    error('Unsupported parallel setting name: %s.',name); end

% look up current setting
try
    res = tracking.parallel.(name);
catch
    defaults = struct('engine',{'BLS'},'pool',{{}},'policy',{'par_reschedule_policy'},'verbosity',{0},'logfiles',{});
    tracking.parallel.(name) = defaults.(name);
    res = tracking.parallel.(name);
end

% optionally override current setting
if nargin >= 2
    % perform sanity checks
    switch name
        case 'engine'
            if ~ischar(val) || ~any(strcmp(val,{'BLS','local','ParallelComputingToolbox','Reference'}))
                error('Unsupported value for the ''engine'' setting: %s.',hlp_tostring(val)); end
        case 'pool'
            if ~iscellstr(val)
                error('The ''pool'' setting must be a cell array of strings.'); end
        case 'policy'
            if ~ischar(val) || ~exist(val,'file')
                error('The ''policy'' setting must be a string an refer to an existing function.'); end
        case 'verbosity'
            if ~isscalar(val) || ~isnumeric(val) || val~=round(val) || val<0
                error('The ''verbosity'' setting must be a nonnegative integer.'); end
        case 'logfiles'
            if ~iscellstr(val)
                error('The ''logfiles'' setting must be a cell array of strings.'); end
        otherwise
            error('Unsupported parallel setting: %s.',name);
    end
    tracking.parallel.(name) = val; 
end
