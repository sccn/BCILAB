function harvested_addresses = par_parse_logfiles(logpaths, harvest_timeout, harvest_ips)
% Parse logfiles produced by par_worker, extract pid@host:port
% Addresses = par_parse_logfiles(LogPaths, Timeout)
%
% This function can scan multiple logfiles simultaneously, wait until the required information has
% been written and extracted, or until a timeout has expired, and then returns the parsed information.
%
% In:
%   LogPaths : cell array of paths to logfiles, or directory name/wildcard pattern, or single
%              file name
%
%   Timeout : timeout in seconds until the function returns with the information found so far
%             (default: 300)
%
%   GrabIPs : whether to grab the IP addresses rather than hostnames (default: false)
%
% Out:
%   Addresses : a cell array of 'pid@host:port' strings, one per worker that was found, in no
%               particular order
%
%                                Christian Kothe, Syntrogi
%                                2015-04-16

% Copyright (C) Christian Kothe, Syntrogi, 2015, christian.kothe@syntrogi.com
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

if ischar(logpaths)
    if exist(logpaths) == 7
        % logpaths is a directory
        records = rdir(logpaths);
        logpaths = {records([records.isdir] == 0).name};
    elseif exist(logpaths) == 2
        logpaths = {logpaths}; 
    else
        error('BCILAB:file_not_found','The file %s was not found.',logpaths);
    end
end
if nargin < 2
    harvest_timeout = 300; end
if nargin < 3
    harvest_ips = false; end

% harvest the pid@host:port information from the log files...
host_line = 'this is bcilab worker';
port_line = 'listening on port';
proc_line = 'the process id is';
error_line = 'no free port found; exiting';
harvested_addresses = {};   % the list of host:port addresses harvested from log files so far
active_logfiles = logpaths; % the log files that are still actively being scanned
t0 = tic;
while toc(t0) < harvest_timeout
    for k=length(active_logfiles):-1:1
        fn = active_logfiles{k};
        if exist(fn,'file')
            try
                fid = fopen(fn,'r');
                if fid==-1
                    continue; end
                content = vec(fread(fid,Inf,'*char'))';
                host_match = strfind(content,host_line);
                port_match = strfind(content,port_line);
                proc_match = strfind(content,proc_line);
                if ~isempty(strfind(content,error_line))
                    % worker had an error: remove it from the set of logfiles being tracked
                    active_logfiles(k) = [];                    
                elseif ~isempty(host_match) && ~isempty(port_match) && ~isempty(proc_match)
                    % found host and port lines
                    host_startofs = host_match(1)+length(host_line); host_endofs = host_startofs + find(content(host_startofs:end)==10,1); host_section = hlp_split(strtrim(content(host_startofs:host_endofs-3)),'/');
                    port_startofs = port_match(1)+length(port_line); port_endofs = port_startofs + find(content(port_startofs:end)==10,1); port_section = strtrim(content(port_startofs:port_endofs-3));
                    proc_startofs = proc_match(1)+length(proc_line); proc_endofs = proc_startofs + find(content(proc_startofs:end)==10,1); proc_section = strtrim(content(proc_startofs:proc_endofs-3));
                    % record harvested address and remove it from the set of tracked logfiles
                    harvested_addresses{end+1} = [proc_section '@' host_section{1+harvest_ips} ':' port_section]; %#ok<AGROW>
                    fprintf('  found %s\n', harvested_addresses{end});
                    active_logfiles(k) = [];
                end
            catch e
                fclose(fid);
                rethrow(e);
            end
        end
    end
    if isempty(active_logfiles)
        break; end
    pause(5);
end

