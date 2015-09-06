function mpw(num)
%   Set the memory per worker that shall be requested from the cluster, in GB.
%
%   Example: wpn 50
%
%                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                               2010-04-22

global tracking;
opts = hlp_varargin2struct(tracking.acquire_options);
if ischar(num)
    num = str2num(num); end
opts.min_memory = num*2^30;
tracking.acquire_options = hlp_struct2varargin(opts);
