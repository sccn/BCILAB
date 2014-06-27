function [signal,state] = set_makepos(varargin)
% Create epochs relative to the target markers in a data set.
% [Signal,State] = set_makepos(Signal,TimeWindow,OnlineEpoching,State)
%
% This function turns a continuous data set with "target" markers (aka target events), i.e. events
% that have a .target field which contains a non-empty "target value" into an epoched data set,
% creating one epoch relative to each target marker, and copying that markers's .target value into the
% respective epoch's .target field. Target markers are defined via the function set_targetmarkers.
%
% In:
%   Signal : continuous data set from which epochs shall be extracted
%
%   TimeWindow : time window relative to the events (forming the extracted epoch), in seconds
%
%   OnlineEpoching : Online epoch extraction scheme; can be one of the following options: 
%                    * 'at_end' : a single epoch is extracted at the end of the signal (default)
%                    * 'at_targets' : one epoch per target marker in the signal is extracted
%
%   State : previous filter state, as obtained by a previous execution of set_makepos on an
%           immediately preceding data set (default: [])
%
% Out:
%   Signal : epoched (segmented) data set that is compatible with EEGLAB's epoch structure:
%            * .data and every other time series field is a 3d array (#channels x #timepoints x #segments)
%            * .xmin/xmax denote the epoch bounds in seconds (relative to epoch-generating event)
%            * .epoch contains the fields of the original epoch-generating events (including
%              .target)
%            * .event contains the flat array of events that overlap the epochs, possibly with
%              replication where .event.latency contains the cumulative latency of all earlier
%              epochs for compatbility with EEGLAB
%            * .epoch.event* for each field .event.* contains a cell array of the field values for
%              all events that overlap the epoch, where .epoch.eventlatency and .epoch.eventduration
%              are in miliseconds rather than samples for compatibility with EEGLAB
%            * .event.epoch contains the index of the epoch that the event lies in and .epoch.event
%              contains an array of indices into .event (all epoch-overlapping events) 
%
%   State :  state of the filter, after it got applied
%
% Notes:
%   In the regular BCILAB processing pipeline, target markers are defined on the continuous 
%   data set either automatically by the bci_ functions, or manually via set_targetmarkers.
%   One can also manually assign a .target field to a subset of markers.
%
% Examples:
%   % extract epochs from a continuous data set, beginning 0.5s before the position of each marker
%   % and ending 1.5s after the position of the corresponding marker
%   eeg = set_makepos(eeg,[-0.5 1.5])
%
%   % as before, but first add target markers for event types 'X' and 'Y'
%   eeg = set_makepos(set_targetmarkers(eeg,{'X','Y'}),[-0.5 1.5])
%
% See also:
%   set_targetmarkers, bci_train
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_makepos_version<2.01> -- for the cache

if ~exp_beginfun('filter'), return; end

declare_properties('name',{'EpochExtraction','epoch'}, 'independent_channels',true,'independent_trials',true);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'time_bounds','TimeWindow','epobounds'}, [], [], 'Epoch time window relative to the target markers. In seconds.','type','expression','shape','row'), ...
    arg_nogui({'online_epoching','OnlineEpoching'}, 'at_end', {'at_end','at_targets'}, 'Online epoch extraction. If set to at_end, an single epoch is extracted at the end of the signal; if set to at_target_markers, one epoch per target marker in the signal is extracted.'), ...
    arg_nogui({'state','State'}));


% handle previous state for incremental processing (if any)
if isempty(state) %#ok<*NODEF>
    % no state yet: determine epoching info
    state.sample_bounds = round(time_bounds*signal.srate);
    state.sample_range = state.sample_bounds(1):state.sample_bounds(2);
    state.sample_points = length(state.sample_range);
    state.time_bounds = state.sample_bounds/signal.srate;    
    % initialize data buffer
    state.timeseries_fields = utl_timeseries_fields(signal);
    for f = state.timeseries_fields
        state.buffer.(f{1}) = []; end
    % initialize marker buffer
    state.events = [];
    % determine how many samples should to be carried over between calls 
    % (actual buffer size may be less if not enough data is available)
    state.prepend_samples = max([state.sample_points-1,state.sample_range(end),-state.sample_range(1)]);
else
    % have state: prepend the buffer to the signal
    for f = state.timeseries_fields
        signal.(f{1}) = cat(2,state.buffer.(f{1}),signal.(f{1})); end
    [signal.nbchan,signal.pnts,signal.trials,extra_dims] = size(signal.data); %#ok<NASGU>
   % increment marker latencies by amount of prepended data
    if ~isempty(signal.event)
        [signal.event.latency] = arraydeal([signal.event.latency]+size(state.buffer.data,2)); end
    % prepend buffered markers
    try
        signal.event = [state.events signal.event];
    catch %#ok<CTCH>
        disp('Event structure changed; discarding buffered events.');
        state.events = [];
    end
end

% update final state
if nargout > 1
    % retain the last portion of the time series fields
    for f = state.timeseries_fields
        state.buffer.(f{1}) = signal.(f{1})(:,max(1,end-state.prepend_samples+1):end,:,:,:,:,:,:); end
    if ~isempty(signal.event)
        % retain events that lie in the buffered interval
        state.events = signal.event(round([signal.event.latency]) > (signal.pnts-size(state.buffer.data,2)));
        % change their latency: the last sample's latency goes from pnts to the state buffer's size
        if ~isempty(state.events)
            [state.events.latency] = arraydeal([state.events.latency] + (size(state.buffer.data,2)-signal.pnts)); end
    else
        % otherwise there are no events to carry over
        state.events = [];
    end
end

if strcmp(online_epoching,'at_targets') || ~onl_isonline
    % default offline case: generating epochs around target markers
    if ~isempty(signal.event)
        % make sure that the necessary fields are present
        if isfield(signal,'epoch') && ~isempty(signal.epoch)
            error('The signal is already epoched; can only extract epochs from continuous data.'); end
        if ~isfield(signal.event,'latency')
            error('The markers of this signal must have a field named .latency (which holds the latency of each marker).'); end
        if ~isfield(signal.event,'target')
            error('The markers of this signal must have a field named .target for epoch generation to work (see also set_targetmarkers).'); end

        % make sure that events are sorted
        latencies = [signal.event.latency];
        if ~issorted(latencies)
            [latencies,order] = sort(latencies);
            signal.event = signal.event(order);
        end
        
        % make sure that all event latencies fall within the signal bounds
        out_of_bounds = round(latencies)<1 | round(latencies)>signal.pnts;
        if any(out_of_bounds)
            signal.event(out_of_bounds) = [];
            latencies(out_of_bounds) = [];
        end
        
        % determine the indices/latencies of events that shall be used for epoch generation
        target_events = find(~cellfun('isempty',{signal.event.target}));
        target_latencies = round(latencies(target_events));
        
        % remove target events whose epoch ranges intersect data boundaries
        out_of_bounds = target_latencies+state.sample_range(1)<1 | target_latencies+state.sample_range(end)>signal.pnts;
        if any(out_of_bounds)
            target_events(out_of_bounds) = []; 
            target_latencies(out_of_bounds) = []; 
        end
        
        if ~isempty(target_events)
            % generate index ranges for epoching
            ranges = bsxfun(@plus,target_latencies,state.sample_range(:));

            % epoch the time series fields
            for f = state.timeseries_fields
                siz = size(signal.(f{1}));
                if ~isempty(signal.(f{1}))
                    signal.(f{1}) = reshape(signal.(f{1})(:,ranges,1,:,:,:,:,:),[siz(1),size(ranges,1),size(ranges,2),siz(4:end)]); end
            end

            % epoch sparse event latencies and collect retained event indices and latencies
            [event_positions,residuals] = sparse_binning(latencies,[],signal.pnts);
            [dummy,sample_indices,event_indices] = find(event_positions(:,ranges(:))); %#ok<ASGLU>
            residuals = residuals(:); sample_indices = sample_indices(:); event_indices = event_indices(:);
            
            % get the associated epoch for each event, and the number of events in each epoch
            epoch_assignment = 1+floor((sample_indices-1)/state.sample_points);
            epoch_numevents = full(sum(logical(sparse_binning(epoch_assignment,[],length(target_events))),1));
            
            % initialize .epoch (from time-locking events) and .event (from retained events)
            signal.epoch = signal.event(target_events);
            signal.event = signal.event(event_indices);
            
            % cross-link epochs and events
            if ~isempty(signal.event)
                [signal.event.epoch] = arraydeal(epoch_assignment); end
            [signal.epoch.event] = chopdeal(1:length(event_indices),epoch_numevents);
            
            % copy all event.* fields to epoch.event*
            for f=fast_setdiff(fieldnames(signal.event)',{'latency','duration','epoch'})
                [signal.epoch.(['event' f{1}])] = chopdeal({signal.event.(f{1})},epoch_numevents); end
            
            % special handling for epoch.eventlatency, and epoch.eventduration (in ms and relative to epoch center)
            epoch_rellatencies = (residuals(event_indices)+sample_indices-1-(epoch_assignment-1)*state.sample_points)+state.sample_range(1);
            [signal.epoch.eventlatency] = chopdeal(num2cell(epoch_rellatencies/signal.srate*1000),epoch_numevents);
            if isfield(signal.event,'duration') 
                durations = {signal.event.duration};
                [durations{cellfun('isempty',durations)}] = deal(0);
                [signal.epoch.eventduration] = chopdeal(num2cell([durations{:}]*(1000/signal.srate)),epoch_numevents);
            end
            
            % update event latency to pseudo-continuous latencies
            if ~isempty(signal.event)
                [signal.event.latency] = arraydeal(sample_indices + residuals(event_indices)); end
        else
            for f = state.timeseries_fields
                signal.(f{1}) = []; end
        end
    else
        for f = state.timeseries_fields
            signal.(f{1}) = []; end        
    end
else
    % default online case: extract a single epoch at the current end of the stream
    if signal.pnts > state.sample_points
        % too long: remove excess data
        for f = state.timeseries_fields
            if ~isempty(signal.(f{1}))
                signal.(f{1}) = signal.(f{1})(:,end-state.sample_points+1:end,1,:,:,:,:,:); end
        end
        % update event latencies and remove out-of-range events
        if ~isempty(signal.event)
            [signal.event.latency] = arraydeal([signal.event.latency]+(state.sample_points-signal.pnts));
            signal.event(round([signal.event.latency])<1) = [];
        end
    elseif signal.pnts < state.sample_points
        % too short: generate no epoch
        for f = state.timeseries_fields
            signal.(f{1}) = []; end
    end
    
    if ~isempty(signal.data)
        % add an epoch field
        signal.epoch = struct('latency',{signal.smax - state.time_bounds(2)*signal.srate});
        if ~isempty(signal.event)
            % cross-link epochs and events
            signal.epoch.event = 1:length(signal.event);            
            [signal.event.epoch] = deal(1);
            % copy event fields over
            for f=fast_setdiff(fieldnames(signal.event)',{'latency','duration','epoch'})
                signal.epoch.(['event' f{1}]) = {signal.event.(f{1})}; end
            % apply special handling to latency and duration (both in ms and relative to epoch center)
            signal.epoch.eventlatency = num2cell(([signal.event.latency]-1+state.sample_range(1))/signal.srate*1000);
            if isfield(signal.event,'duration')
                durations = {signal.event.duration};
                [durations{cellfun('isempty',durations)}] = deal(0);                
                signal.epoch.eventduration = num2cell([durations{:}]*(1000/signal.srate)); 
            end
        end
    end
end

% fix up signal properties
signal.xmin = state.time_bounds(1);
signal.xmax = state.time_bounds(2);
[signal.nbchan,signal.pnts,signal.trials,extra_dims] = size(signal.data); %#ok<NASGU>
if isempty(signal.data) || isempty(signal.event)
    signal.event = []; end
if isempty(signal.epoch)
    signal.epoch = struct([]); end

exp_endfun;
