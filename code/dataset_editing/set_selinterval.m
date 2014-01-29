function signal = set_selinterval(varargin)
% Selects a time interval from a data set.
% Signal = set_selinterval(Signal,Intervals,Unit)
%
% In:
%   Signal    :   Continuous EEGLAB data set from which an interval should be selected
%
%   Intervals :   Selection interval(s) formatted as [start end; start end; start end; ...] or
%                 a MATLAB index range for the samples (if IntervalUnit is set to 'range')
%
%   IntervalUnit : Unit of measurement for the interval. Either 'seconds', 'samples', 'fraction',
%                  or 'range' (default: seconds)
%
%   InsertBoundaryMarkers : Whether to insert boundary markers (default: false)
%
% Out:
%   Signal    :   data set restricted to the selected range. The following fields are updated:
%                 * .data and all other time series fields are reduced to the selected intervals
%                 * .event is reduced to selected events
%                 * .xmax/.pnts are updated
%                 * optionally boundary events are inserted
%
% Examples:
%   % for a continuous data set, retain only the data within 50s and 200s, as well as 1200s and 1500s
%   eeg = set_selinterval(eeg,[50 200; 1200 1500])
%
%   % for a continuous data set, retain only the data within the 1000's and the 10000's sample
%   eeg = set_selinterval(eeg,[1000 10000],'samples')
%
%   % for a continuous data set, retain only the last half of the data
%   eeg = set_selinterval(eeg,[0.5 1],'fraction')
%
%   % as before, but pass the arguments by name
%   eeg = set_selinterval('Signal',eeg, 'Intervals',[0.5 1], 'Unit','fraction')
%
% See also:
%   set_selepos, set_partition
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01

% set_selinterval_version<2.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'intervals','Intervals'},[],[], 'Selection intervals. Array of the form [start end; start end; start end; ...].'), ...
    arg({'selunit','IntervalUnit','Unit','unit'},'seconds',{'seconds','samples','fraction','range'}, 'Interval unit. The unit of measurement for selection intervals (does not apply to SampleRange).'), ...
    arg({'insert_boundary_markers','InsertBoundaryMarkers'},false,[], 'Insert boundary markers. Whether to insert boundary markers (for EEGLAB compatibility).'));

if ~isempty(signal.epoch) || size(signal.data,3)>1
    error('This function cannot be applied to epoched data. Use flt_window to apply window functions to epoched data.'); end

switch selunit
    case 'fraction'
        intervals = 1+round(intervals*(size(signal.data,2)-1));
        samplerange = [];
    case 'seconds'
        intervals = min(1+round(intervals*signal.srate),size(signal.data,2));
        samplerange = [];
    case 'range'
        samplerange = intervals(:)';
        intervals = [];
    case 'samples'
        samplerange = [];
    otherwise
        error(['Unsupported selection unit: ' selunit]);
end

% generate samplerange from intervals if not yet present
if isempty(samplerange)
    samplerange = false(1,size(signal.data,2));
    for k=1:size(intervals,1)
        samplerange(intervals(k,1):intervals(k,2)) = true; end
end

if ~any(samplerange)
    % clear data
    for field = utl_timeseries_fields(signal)
        signal.(field{1}) = []; end
    signal.event = [];
    signal.pnts = 0;
elseif length(samplerange)<size(signal.data,2) || ~all(samplerange)
    % select range within the time series fields
    for field = utl_timeseries_fields(signal)
        if ~isempty(signal.(field{1}))
            signal.(field{1}) = signal.(field{1})(:,samplerange,:,:,:,:,:,:); end
    end
    % select range within the events
    if ~isempty(signal.event)
        % bin latencies into sparse array and apply selection
        [event_positions,residuals] = sparse_binning([signal.event.latency],[],signal.pnts);
        [ranks,sample_indices,event_indices] = find(event_positions(:,samplerange));
        residuals = residuals(:); sample_indices = sample_indices(:); event_indices = event_indices(:);
        % write event subset and update latencies
        if ~isempty(ranks)
            signal.event = signal.event(event_indices);
            [signal.event.latency] = arraydeal(sample_indices+residuals(event_indices)); 
        else
            signal.event = [];
        end
    end
    % optionally insert boundary markers
    if insert_boundary_markers
        % generate data intervals
        if isempty(intervals)
            sample_mask = logical([]);
            sample_mask(samplerange) = true;
            intervals = reshape(find(diff([false sample_mask false])),2,[])';
            intervals(:,2) = intervals(:,2)-1;
        end
        % append new events for each boundary between intervals
        range = length(signal.event)+(1:size(intervals,1)-1);
        if range
            [signal.event(range).type] = deal('boundary');
            [signal.event(range).latency] = arraydeal(intervals(1:end-1,2));
            [signal.event(range).duration] = arraydeal(intervals(2:end,1) - (intervals(1:end-1,2)+1));                    
            % resort events
            [dummy,order] = sort([signal.event.latency]); %#ok<ASGLU>
            signal.event = signal.event(order);
        end
    end
    % update misc fields
    signal.pnts = size(signal.data,2);
    signal.xmax = signal.xmin + (find(samplerange,1,'last')-1)/signal.srate;
    signal.xmin = signal.xmin + (find(samplerange,1)-1)/signal.srate;
end

exp_endfun;
