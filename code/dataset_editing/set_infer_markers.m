function signal = set_infer_markers(varargin)
% Infers markers from an event channel, if possible.
% Signal = set_infer_markers(Signal)
%
% If the data set has one (or more) channels with relatively few unique values, and relatively few
% off-median values (or plateaus thereof), markers are being generated to encode whenever such a
% channel changes its value. By default, the detected event channels are being removed from the
% data. 
%
% This funciton is automatically called by io_loadset if no markers are present in the given 
% set. io_loadset has a parameter ('markerchannel') which can be utilized to customize the settings 
% below (as a cell array of arguments to set_infer_markers).
%
% In:
%   Signal                    : continuous EEGLAB data set from which an interval should be selected
%
%   MaxEvents                 : a channel is only considered for event generation if it would produce
%                               at most this many distinct events (default: 25000)
%
%   MaxTypes                  : a channel is only considered for event generation if it would produce
%                               at most this many distinct event types (default 300)
%
%   IncludeLabel              : when to include the channel label into the generated event types; 
%                               can be 'always', 'multiplechans', 'never' (default: 'multiplechans')
%
%   AllPositive               : whether all events must have positive offsets from the baseline 
%                               (median) value (default: true)
%
%   EncodePlateaus            : Whether the trigger channel may contain (non-baseline) plateaus, to 
%                               be encoded in event duration. (default: true)
%
%   StrictInteger             : only consider channels that would produce integer event types 
%                               (otherwise: round to integer) (default: true)
%
%   RelativeBaselineThreshold : Relative baseline threshold. If the event types would be larger than
%                               this threshold, they will be re-coded to be relative to the smallest
%                               event type in the channel. (default: 10000)
%
%   OmitBaselineThreshold     : Baseline omission threshold. If the baseline event takes up more than
%                               this fraction of the data, it will not be encoded as events.
%                               (default: 0.8)
%
%   RemoveEventChannels       : Whether to remove channels that have been identified as event
%                               channels (default: true)
%
% Out:
%   Signal    : data set restricted to the selected range
%
% Examples:
%   % for a data set with one or more trigger channels, fill in contents of the .event field
%   % and remove the trigger channel(s); do nothing if the data set does not contain trigger channels
%   eeg = set_infer_markers(eeg)
%
%   % fill in contents of the .event field from trigger channels, and do not remove these channels
%   eeg = set_infer_markers('Signal',eeg,'RemoveEventChannels',false)
%
%   % if an event channel is not being detected because it produces more than the default MaxEvents 
%   % number of events, use a larger cutoff value
%   eeg = set_infer_markers('Signal',eeg,'MaxEvents',30000)
%   
%   % if an event channel is not being detected because it produces more than the default MaxTypes
%   % number of distinct event types, use a larger cutoff value
%   eeg = set_infer_markers('Signal',eeg,'MaxTypes',1000)
%
%   % for event channels that have values that far from zero (say, 30000+), the values are by default
%   % made relative to to the baseline value (most frequent or smallest value, depending on the 
%   % EncodePlateaus setting) -- if this is not intended, override the cutoff with a larger value
%   eeg = set_infer_markers('Signal',eeg,'RelativeBaselineThreshold',Inf)
%
%   % if a trigger channel contains many subsequent events of the same type (i.e. plateaus), by 
%   % default only one event is generated with the appropriate duration; if instead an event should 
%   % be generated for each sample, override that default (and possibly increase MaxEvents)
%   eeg = set_infer_markers('Signal',eeg,'EncodePlateaus',false,'MaxEvents',Inf)
%   
% See also:
%   io_loadset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-01-26

declare_properties('independent_channels',false,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'force_processing','ForceProcessing'},false,[], 'Force marker processing. Scan for marker channels even if a data set already has events.'), ...
    arg({'max_events','MaxEvents'},30000,uint32([0 1000000000]), 'Number of events allowed. If a channel would produce more than this many events, it is not considered an event channel.'), ...
    arg({'max_event_fraction','MaxEventFraction'},0.3,[0 1], 'Maximum fraction of events. If a channel would produce more events than this fraction of its total number of samples, the channel is not considered an event channel.'), ...
    arg({'max_types','MaxTypes'},300,uint32([0 1000000000]), 'Number of eventtypes allowed. If a channel would produce more than this many event types, it is not considered an event channel.'), ...
    arg({'include_label','IncludeLabel'},'multiplechans',{'always','multiplechans','never'},'Integrate channel labels. When to integrate channel labels into the event types - by default. only when there are multiple event channels.'), ...
    arg({'all_positive','AllPositive'},true,[],'Events must be positive. Only if all event types produced by a channel would correspond to positive numbers, the channel is considered an event channel.'), ...
    arg({'encode_plateaus','EncodePlateaus'},true,[],'Plateaus encoded. The trigger channel may contain (non-baseline) plateaus, to be encoded in event duration.'), ...
    arg({'strict_integer','StrictInteger'},true,[],'Strictly integer. Only allow channels that have integer event types (otherwise: round to integers).'), ...
    arg({'relative_baseline_thresh','RelativeBaselineThreshold'},10000,[],'Relative baseline threshold. If the event types would be larger than this threshold, they will be re-coded to be relative to the smallest event type in the channel.'), ...
    arg({'omit_baseline_thresh','OmitBaselineThreshold'},0.8,[0 1],'Baseline omission threshold. If the baseline event takes up more than this fraction of the data, it will not be encoded as events.'), ...
    arg({'remove_eventchns','RemoveEventChannels'},true,[],'Remove event channels. Whether to remove those channels that have been identified as event channels.'));

% obtain the chanlocs
if isfield(signal,{'head','parts'})
    signal = exp_eval(signal); end
utl_check_fields(signal,{'data','chanlocs'},'signal','signal');

% figure out which are the event channels
fprintf('Scanning potential marker channel ');
for k=1:size(signal.data,1) %#ok<*NODEF>
    fprintf('%i ',k);
    X = signal.data(k,:);
    X(isnan(X)) = 0;
    % do a few vectorized computations for speed...
    numtypes(k) = length(unique(X));
    if numtypes(k) <= max_types
        if encode_plateaus
            % in the plateaus case, the baseline is 0 or the smallest value
            baseline(k) = min(X);
            if baseline(k) < relative_baseline_thresh && baseline(k) > 0
                baseline(k) = 0; end
            % whether baseline periods are encoded as events depends on the fraction of time in the
            % data set
            omit_baseline(k) = mean(X == baseline(k)) > omit_baseline_thresh;
            eventmask = ([1; diff(X(:))] ~= 0);
            if omit_baseline(k)
                eventmask = eventmask & (X(:) ~= baseline(k)); end
            numevents(k) = nnz(eventmask);
        else
            % without plateaus, the baseline type is the mode
            baseline(k) = mode(X);
            numevents(k) = nnz(X - baseline(k));
            omit_baseline(k) = true;
        end
    else
        % channel will not be considered
        baseline(k) = 0;
        omit_baseline(k) = true;
        numevents(k) = 0;
    end
    allpositive(k) = all(X >= baseline(k));
    allinteger(k) = all(abs(X-round(X)) <= eps(X));
end
fprintf('\n');

% find the mask of event channels
eventchans = (numtypes <= max_types) & (numevents <= max_events) & (numevents <= size(signal.data,2)*max_event_fraction) & (~all_positive | allpositive) & (~strict_integer | allinteger);
    
% determine whether to include labels
include_label = hlp_rewrite(include_label,'always',true,'never',false,'multiplechans',nnz(eventchans) > 1);

% initialize events
if ~isfield(signal,'event') || isempty(signal.event)
    signal.event = struct('latency',{},'duration',{},'type',{}); end

% generate events
for k = find(eventchans)
    try
        % find event codes per sample
        X = round(signal.data(k,:) - baseline(k));
        X(isnan(X)) = 0;
        % do a run-length encode into events
        lat = find([1; diff(X(:))] ~= 0);
        % optionally mask out baseline events
        if omit_baseline(k)
            mask = X(lat) ~= 0;
        else
            mask = true(1,length(lat));
        end
        % get values, duration, and latency as cell arrays
        dur = num2cell(diff([lat; numel(X)+1]));
        val = cellfun(@num2str,num2cell(X(lat(mask)),1),'UniformOutput',false);
        if include_label
            val = cellfun(@(n) [signal.chanlocs(k).labels '_' n],val,'UniformOutput',false); end
        lat = num2cell(lat(mask));
        dur = dur(mask);
        if ~isempty(lat)
            % make space
            signal.event(end+length(lat)).latency = [];
            % insert new content
            [signal.event(end-length(lat)+1:end).latency] = lat{:};
            [signal.event(end-length(lat)+1:end).duration] = dur{:};
            [signal.event(end-length(lat)+1:end).type] = val{:};
        else
            % the channel is not actually an event channel (but a flatline)
            eventchans(k) = false;
        end
    catch e
        disp(['Could not process channel ' num2str(k)]);
        env_handleerror(e);
    end
end

% sort events
signal.event = signal.event(hlp_getresult(2,@sort,[signal.event.latency]));

% remove event channels
signal.data = signal.data(~eventchans,:,:,:,:,:,:,:);
signal.chanlocs = signal.chanlocs(~eventchans);
signal.nbchan = size(signal.data,1);
