function wpn(num)
%   Set the workers per node that shall be requested from the cluster.
%
%   Example: wpn 4
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22

global tracking;
opts = hlp_varargin2struct(tracking.acquire_options);
if isfield(opts,'processors_per_node')
    if ischar(num)
        num = str2num(num); end
    opts.processors_per_node = num;
elseif isfield(opts,'num_workers')
    fprintf('Only num_workers can be set in current acquire mode. Setting that instead.\n')
    opts.num_workers = num;
end
tracking.acquire_options = hlp_struct2varargin(opts);
