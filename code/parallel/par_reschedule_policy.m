function result = par_reschedule_policy(batchid,inflight,waiting,times,msgs)
% Default re-scheduling policy for the Java-based Scheduler.
% Reschedule = par_reschedule_policy(Batch-Id,Inflight-Ids,Event-Times,Event-Messages)
%
% This function is periodically invoked by the Scheduler with a list of recent events (times in ms
% and content encoded as strings) and tags of in-flight tasks. It is expected to issue the
% re-scheduling of starved or lost tasks (out of those that are in-flight), depending on some
% assumptions.
%
% In:
%   Batch-Id     : Identifier of the current batch of tasks (this policy may be invoked for multiple
%                  possibly overlapping schedules).
%
%   Inflight-Ids : Cell array of Ids/Tags of in-flight tasks (same as in the Event-Msgs)
%
%   Waiting-Ids : Cell array of Ids/Tags of waiting tasks (same as in the Event-Msgs)
%
%   Event-Times  : Cell array of event timestamps in miliseconds (since beginning of the Scheduler's
%                  current batch of tasks).
%
%   Event-Messages: Cell array of event content, indexed like Event-Times.
%
% Out:
%   Reschedule : Java Vector of Integers, referring to the Inflight-Ids that shall be rescheduled
%                (or the empty Vector).
%
% Notes:
%   A log of what was rescheduled and when is collected in the global variable 
%   tracking.temp.rescheduled.
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

global tracking;

% set this to true if your cluster is reliable and you want to minimize resources
no_reschedules = false;

if no_reschedules
    import java.lang.*;
    import java.util.*;
    result = java.util.Vector();
    return;
end

try
        
    % here we keep track of our meta-data for each of the batches
    persistent batches; %#ok<TLEV>
    default_batch = struct('rescheduled',{[]});  % set of already rescheduled tasks for this batch
    if length(batches) < batchid
        if isempty(batches)
            batches = default_batch; end
        batches(batchid) = default_batch;
    end
    
    % turn times into seconds relative to the first event
    times = [times{:}]/1000;
    starttime = min(times);
    times = (times-starttime);
    
    % turn inflight into a vector
    inflight = [inflight{:}];
    
    
    
    % --- begin policies ---
    reschedule = [];
    
    
    % Policy 1: endgame mode: if less than N tasks in flight, re-schedule them to other machines (so
    % that they are being worked on by two machines)
    if isempty(waiting) && length(inflight) <= 1
        % but do not reschedule them twice
        reschedule = setdiff(inflight,batches(batchid).rescheduled);
    end
    
    
    % Policy 2: identify stragglers, and re-schedule them, assuming that task completion time is
    % normally distributed
    try
        % get a table of recorded tasks
        tasks = struct();
        for m=1:length(msgs)
            msg = msgs{m};
            assigned = strncmp('assigned',msg,8);
            received = strncmp('received',msg,8);
            if assigned || received
                tmp = hlp_split(msg,':');
                taskid = ['x' tmp{2}];
                if assigned
                    tasks.(taskid).assigned = times(m);
                elseif received
                    tasks.(taskid).received = times(m);
                end
                tasks.(taskid).worker = tmp{3};
            end
        end
        
        % get the task completion times for the completed tasks
        completion_times = [];
        for f=fieldnames(tasks)'
            t = tasks.(f{1});
            if isfield(t,{'assigned','received'})
                completion_times(end+1) = t.received - t.assigned; end
        end
        
        % if estimates are reasonable at this point (half of the jobs have been scheduled)
        if length(completion_times) >= length(inflight) && length(completion_times) >= 5
            tnow = double(tic)/1000000 - starttime;
            % estimate the parameters mu,sigma assuming a truncated normal distribution of
            % completion times
            custompdf = @(x,mu,sigma) normpdf(x,mu,sigma)./normcdf(tnow,mu,sigma);
            
            [estim,conf] = mle(completion_times,'pdf',custompdf,'start',[mean(completion_times),std(completion_times)],'lower',[-Inf 0]);       %#ok<NASGU>
            mu = estim(1);
            sigma = estim(2);
            
            % for all inflight (and registered) tasks...
            for i=1:length(inflight)
                t = inflight(i);
                taskid = ['x' num2str(t)];
                % estimate their duration (tasks for which we have no record are assumed to have
                % started before our records began)
                if isfield(tasks,taskid) && isfield(tasks.(taskid),'asssigned')
                    duration(i) = tnow - tasks.(taskid).assigned;
                else
                    duration(i) = tnow;
                end
            end
            
            % decide which ones we reschedule, based on how long other tasks took so far
            for i=1:length(inflight)
                taskid = inflight(i);
                % duration relative to mean, in std. devs (mahalanobis metric)
                mahal = (duration(i) - mu) / sigma;
                if mahal > 3 && isempty(find(batches(batchid).rescheduled,taskid))
                    % is at least >3 std devs, and the task has not yet been rescheduled
                    reschedule(end+1) = taskid;
                end
            end
        end
        
    catch e
        if ~any(strcmp(e.identifier,{'stats:mle:NonfinitePdfVal','stats:mle:NonpositivePdfVal'}))
            disp('Issue in the reschedule policy: ');
            hlp_handleerror(e);
        end
    end
    
    
    % --- end policies ---
    
    % log what is being rescheduled and when
    if ~isempty(reschedule)
        try
            tracking.temp.rescheduled{end+1} = [now reschedule];
        catch
            tracking.temp.rescheduled = {[now reschedule]};
        end
    end
    
    reschedule = unique(reschedule);
    if ~isempty(reschedule)
        batches(batchid).rescheduled = [batches(batchid).rescheduled reschedule]; end
    
    % emit the result as a Java Vector (pre-Generics)
    import java.lang.*;
    import java.util.*;
    result = java.util.Vector();
    if ~isempty(reschedule)
        for k=reschedule
            result.add(Integer(k)); end
    end
    
catch e
    disp('Error executing the reschedule policy; traceback:');
    hlp_handleerror(e);
    result = java.util.Vector();
end