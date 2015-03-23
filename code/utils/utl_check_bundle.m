function bundle = utl_check_bundle(bundle)
% Check a stream bundle for consistency and fix if necessary.
%
% This function ensures that:
%  * all streams cover the same time interval (up to respective sample accuracy)
%  * each stream has all markers (with matching type, latency and target field; each stream may have 
%    its unique additional marker meta-data)
%  * the expression associated with each stream (if any) matches the data (via a stored fingerprint)
%    unless that functionality is disabled
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
dp;

% turn into a bundle if necessary
if ~isfield(bundle,'streams')
    bundle = struct('streams',{{bundle}}); end

str = bundle.streams;

% some aggregate information
mrk_type = {};
mrk_latency = {};
mrk_target = {};
mrk_source = [];
xmin = [];
xmax = [];

% evaluate and validate inputs
for s=1:length(str)
    % evaluate if necessary
    str{s} = exp_eval_optimized(str{s});
    utl_check_fields(str{s},{'data','epoch','event','xmin','xmax','srate','pnts','nbchan','trials','chanlocs'},'signal','signal');
    % validate event field
    if ~isempty(str{s}.event)
        if ~isfield(str{s}.event,'type')
            error('WARNING: One of the signals is lacking the .event.type field. This likely leads to errors during later processing.'); end
        if ~all(cellfun('isclass',{str{s}.event.type},'char') | cellfun('isempty',{str{s}.event.type}))
            error('The signal contains events whose .type field is non-empty and not a string. This is an error.'); end
        if ~isfield(str{s}.event,'latency')
            error('One of the signals is lacking the .event.latency field.'); end
        if ~isempty(str{s}.epoch) && ~isfield(str{s}.event,'epoch')
            error('One of the signals is lacking the .event.epoch field.'); end
        latencies = {str{s}.event.latency};
        latency_numels = cellfun('prodofsize',latencies);
        if any(latency_numels == 0)
            error('One or more of the events in the given signal have an empty .latency field, which is not permitted.'); end
        if any(latency_numels ~= 1)
            error('One or more of the events in the given signal have a .latency value that is not a scalar, which is not permitted.'); end
        latencies = [latencies{:}];
        out_of_bounds = latencies < 1 | latencies > size(str{s}.data,2);
        if any(out_of_bounds)
            disp_once('WARNING: The signal contains %i events with out-of-bounds latencies.',nnz(out_of_bounds)); end
        if ~isequal(sort(latencies),latencies)
            disp_once('Warning: The signal has unsorted event latencies.'); end
    end
    % validate epoch field
    if isempty(str{s}.epoch)
        if size(str{s}.data,3) > 1
            error('One of the signals has an epoched .data field but an empty .epoch field.'); end
    else
        if length(str{s}.epoch) ~= size(str{s}.data,3)
            error('The number of epochs in the signal (%i) does not match the length of the .epoch field (%i).',length(str{s}.epoch),size(str{s}.data,3)); end
        if ~isfield(str{s}.epoch,'event')
            error('One of the signals in the bundle is lacking the .epoch.event field.'); end
        eventepochs = {str{s}.event.epoch};
        epoch_numels = cellfun('prodofsize',eventepochs);
        if any(epoch_numels == 0)
            error('One or more of the events in the given signal have an empty .epoch field, which is not permitted.'); end
        if any(epoch_numels ~= 1)
            error('One or more of the events in the given signal have an .epoch value that is not a scalar, which is not permitted.'); end
        eventepochs = [eventepochs{:}];
        if any(eventepochs < 1 | eventepochs > size(signal.data,3))
            error('The given signal has invalid (out-of-bounds) .event.epoch indices.'); end
        epochevents = [str{s}.epoch.event];
        if any(epochevents < 1 | epochevents > length(str{s}.event))
            error('The given signal has invalid (out-of-bounds) .epoch.event indices.'); end
    end
    
    % check if the bundle contains a non-trivial online expression
    if isfield(str{s},'tracking') && isfield(str{s}.tracking,'online_expression')
        if ~strcmp(char(str{s}.tracking.online_expression.head),'rawdata')
            % in this case we forget about all processing applied to the dataset and treat it as if it were "raw"            
            disp_once('WARNING: The given dataset has non-trivial BCILAB filters applied to it. Such filters should be applied in the approach instead, and will not be reflected in the model:\n         %s',hlp_tostring(str{s}.tracking.expression,1000));
            str{s} = rmfield(str{s},'tracking');
        end
    end
    
    % some more sanity checks
    if size(str{s}.data,1) ~= length(str{s}.chanlocs)
        error('The number of channels in the .data field (%i) does not match the length of the .chanlocs field (%i).',size(str{s}.data,1),length(str{s}.chanlocs)); end
    if size(str{s}.data,1) ~= str{s}.nbchan
        disp_once('WARNING: The number of channels in the signal''s .data field (%i) does not match the .nbchan value.',size(str{s}.data,1),str{s}.nbchan); end
    if size(str{s}.data,2) ~= str{s}.pnts
        disp_once('WARNING: The number of time points in the signal''s .data field (%i) does not match the .pnts value.',size(str{s}.data,2),str{s}.pnts); end
    if size(str{s}.data,3) ~= str{s}.trials
        disp_once('WARNING: The number of epochs in the signal''s .data field (%i) does not match the .trials value.',size(str{s}.data,3),str{s}.trials); end
    if str{s}.xmax ~= str{s}.xmin + (str{s}.pnts-1)/str{s}.srate
        disp_once('WARNING: The given signal''s .xmin (%.2f), .xmax (%.2f), .pnts (%i) and .srate (%.1f) values are not mutually consistent.',str{s}.xmin,str{s}.xmax,str{s}.pnts,str{s}.srate); end
    if ischar(str{s}.data)
        error('The given signal has non-numeric data (apparently a reference to a file: %s, make sure that you import it correctly)',str{s}.data); end
    if ~isa(str{s}.data,'double')
        disp_once('WARNING: The given data is not in double-precision format; this can cause severe numerical errors.'); end
end

% find intervals that are covered by all streams
for s=length(str):-1:1
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
        str{s} = exp_eval(set_selinterval(str{s},[xmin xmax],'seconds')); end;
end

% collect markers across all streams
for s=length(str):-1:1
    % collect marker types, latencies and target values across streams
    if ~isempty(str{s}.event)
        mrk_type = [mrk_type {str{s}.event.type}]; %#ok<AGROW>
        mrk_latency = [mrk_latency {str{s}.event.latency}];     %#ok<AGROW>
        if ~isfield(str{s}.event,'target')
            str{s}.event(1).target = []; end
        mrk_target = [mrk_target {str{s}.event.target}]; %#ok<AGROW>
        mrk_source = [mrk_source s*ones(1,length(str{s}.event))]; %#ok<AGROW>
    end
end

% map all markers onto their unique identifier (latency_type)
uids = cellfun(@(a,b)sprintf('%i_%s',b,a),mrk_type,mrk_latency,'UniformOutput',false);
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
    str{s} = exp_eval(utl_check_dataset(str{s}));
end

% finalize
bundle.streams = str;

