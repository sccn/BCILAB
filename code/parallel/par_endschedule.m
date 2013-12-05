function [results,errors] = par_endschedule(sched,varargin)
% Wait for completion of a scheduling operation and return results and errors.
% [Results, Errors] = par_endschedule(Id, Options...)
%
% In:
%   Id : scheduler id, obtained from par_beginschedule
%
%   Options... : optional name-value pairs;
%               'keep': keep this scheduler alive for later re-use (default: false)
%                       if false, the scheduler will be destroyed after use, and re-created during the next run
%
%               'spin_interval' : the period at which par_endschedule will check if the scheduler
%                                 has finished its job (default: 0.1)
%
% Out:
%   Results : cell array of results of the scheduled computations
%   Errors  : cell array of {position,exception struct} for those results that could not be evaluated
%             the position may also be unknown in case of more severe errors
%
% See also:
%   par_beginschedule, par_worker, par_schedule
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

results = {};
errors = {};

% read options
opts = hlp_varargin2struct(varargin,'keep',false, 'spin_interval',0.25);

if iscell(sched)
    % locally computed results
    results = cell(1,length(sched));
    for r=1:length(sched)
        results{r} = sched{r}{2}; end
else
    if isfield(sched,'ReferenceResults')
        % collect results from the Reference implementation
        raw = sched.ReferenceResults;
    else
        % collect results from the BLS scheduler
        
        % wait for the scheduler to finish (note: we cannot wait on a condition variable here,
        % as we need the MATLAB thread to be active for managing the reschedule policy)
        while (~sched.sched.done())
            pause(opts.spin_interval); end
        
        % obtain raw results & convert to cell-string array
        strings = sched.sched.results();
        raw = cell(1,length(strings));
        for k=1:length(strings)
            raw{k} = char(strings(k)); end
        
        % terminate scheduler
        if opts.keep
            sched.sched.clear();
        else
            sched.sched.terminate();
        end
    end
    
    % deserialize & reorder the string-formatted results
    for r=1:length(raw)
        try
            % deserialize result
            raw{r} = hlp_deserialize(fast_decode(raw{r}));
            if isfield(raw{r}{2},{'message','identifier','stack'})
                % append to errors
                errors{end+1} = raw{r};
            else
                % put into results
                results{raw{r}{1}} = raw{r}{2};
            end
        catch e
            errors{end+1} = {NaN,e};
        end
    end    
end

% throw the first error, if not requested as cell array
if ~isempty(errors) && nargout <= 1
    rethrow(errors{1}{2}); end
