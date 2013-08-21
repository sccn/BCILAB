function bundle = utl_check_bundle(bundle)
% Check a stream bundle for consistency and fix if necessary.
%
% This function ensures that:
%  * all streams cover the same time interval (up to respective sample accuracy)
%  * each stream has all markers (with matching type, latency and target field; each stream may have 
%    its unique additional marker meta-data)
%  * the expression associated with each stream (if any) matches the data (via a stored fingerprint)
% 
% In:
%   Bundle : stream bundle to be checked; 
%            if this is an EEGLAB data set, it will be turned into a bundle
%
% Out:
%   Bundle : adjusted bundle
%
% See also:
%   utl_check_dataset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-08-28

% turn into a bundle if necessary
if ~isfield(bundle,'streams')
    bundle = struct('streams',{{bundle}}); end

str = bundle.streams;

mrk_type = {};
mrk_latency = {};
mrk_target = {};
mrk_source = [];
xmin = [];
xmax = [];

% find intervals that are covered by all streams
for s=1:length(str)
    % evaluate if necessary
    str{s} = exp_eval_optimized(str{s});
    % check if epoched
    epoched(s) = isempty(str{s}.epoch);
    % get xmin / xmax
    xmin(s) = str{s}.xmin;
    xmax(s) = str{s}.xmax;
end

if length(unique(epoched)) > 1
    error('The streams in a bundle must either all be epoched or all be non-epoched.'); end

% restrict to the intervals that are covered by all streams
xmin = max(xmin);
xmax = min(xmax);
for s=1:length(str)    
    if xmin > str{s}.xmin || xmax < str{s}.xmax
        str{s} = pop_select(str{s},'time',[xmin xmax]); end
end

% collect markers across all streams
for s=1:length(str)
    % collect marker types, latencies and target values across streams
    if ~isempty(str{s}.event)
        mrk_type = [mrk_type {str{s}.event.type}];
        mrk_latency = [mrk_latency {str{s}.event.latency}];    
        if ~isfield(str{s}.event,'target')
            str{s}.event(1).target = []; end
        mrk_target = [mrk_target {str{s}.event.target}];
        mrk_source = [mrk_source s*ones(1,length(str{s}.event))];
    end
end

% map all markers onto their unique identifier (latency_type)
uids = cellfun(@(a,b)[num2str(b) '_' a],mrk_type,mrk_latency,'UniformOutput',false);
% for each stream..
for s=1:length(str)
    % get the list of markers indices that need to be added to it
    [dummy,add_idx] = setdiff(uids,uids(mrk_source==s)); %#ok<ASGLU>
    % append those markers to the stream
    evt = [str{s}.event];
    range = length(evt)+1:length(evt)+length(add_idx);
    [evt(range).type] = mrk_type{add_idx};
    [evt(range).latency] = mrk_latency{add_idx};
    [evt(range).target] = mrk_target{add_idx};
    % re-sort the events by latency
    [dummy,sort_idx] = sort([evt.latency]); %#ok<ASGLU>
    str{s}.event = evt(sort_idx);
    % check for internal consistency of each data set
    str{s} = utl_check_dataset(str{s});
end

% finalize
bundle.streams = str;