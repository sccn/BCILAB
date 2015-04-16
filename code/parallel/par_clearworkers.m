function par_clearworkers(endpoints,nokillall)
% Clear workers from the listed hostnames (or if those are unspecified, from the current pool).
% par_clearworkers(Hostnames,NoKillall)
%
% Note: this function erases all of the user's MATLAB processes on the given nodes; if there are
% unrelated MATLABs running, these will also be killed. This may also shoot down the master process
% itself if it is located on one of the unfortunate machines...
%
% In:
%   Hostnames : cell array of hostnames or 'pid@host' or 'pid@host:port' strings on which to kill
%               workers
%
%   NoKillall : if given, this function will not perform "killall" for hosts where no process id
%               is given (default: false)
%
% Example:
%   % clear all workers on the given machines (of this user)
%   par_clearworkers({'computer1','computer2','computer3'})
%
%   % clear only specific workers with given pid's, and all MATLAB instances on computer3
%   par_clearworkers({'31443@computer1','41443@computer2','computer3'})
%
%   % same as before, but ignore computer3 (no killall)
%   par_clearworkers({'31443@computer1','41443@computer2','computer3'},true)
%
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

% optionally look up endpoints from the global pool
if ~exist('endpoints','var') || isempty(endpoints)
    endpoints = par_globalsetting('pool'); end
if ~iscellstr(endpoints)
    error('The given Hostnames argument must be a cell array of strings.'); end
if ~exist('nokillall','var') || isempty(nokillall)
    nokillall = false; end

% drop any port assignments
for i=1:length(endpoints)
    colons = endpoints{i}==':';
    if any(colons)
        hostpid{i} = endpoints{i}(1:find(colons)-1); 
    else
        hostpid{i} = endpoints{i};
    end
end

% also split off the pid parts
for i=1:length(hostpid)
    match = hostpid{i}=='@';
    if any(match)
        hosts{i} = hostpid{i}(find(match,1)+1:end); 
        pids{i} = hostpid{i}(1:find(match,1)-1);    
    else
        hosts{i} = hostpid{i};
        pids{i} = 'all';
    end
end

% remove duplicates
uniquehosts = unique(hosts);
pids_per_host = repmat({{}},1,length(uniquehosts));
killall = false(1,length(uniquehosts));

% for each unique host, accumulate the pids to kill (or record 'all' if no pid given)
for i=1:length(hosts)    
    % position of this host in uniquehosts, pids_per_host, and killall
    pos = find(strcmp(uniquehosts,hosts{i}));    
    if strcmp(pids{i},'all')
        killall(pos) = true;
    else
        pids_per_host{pos}{end+1} = pids{i};
    end
end

% now clear them
fprintf('Clearing workers on %i machines...',length(uniquehosts));
for h = 1:length(uniquehosts)
    if ~isempty(pids_per_host{h})
        [status,info] = system(sprintf('ssh -x %s kill %s',uniquehosts{h},sprintf('%s ',pids_per_host{h}{:}))); end
    if ~nokillall && killall(h)
        [status,info] = system(sprintf('ssh -x %s killall MATLAB',uniquehosts{h})); end
end %#ok<NASGU,ASGLU>
fprintf('done.\n');
