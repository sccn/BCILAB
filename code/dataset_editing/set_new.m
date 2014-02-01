function signal = set_new(varargin)
% Create a new EEGLAB data set from individual fields.
% Dataset = set_new(Arguments...)
%
% In:
%   Fields  : Pairs of field names and field values to add to the data set. fields not specified are
%             taken from eeg_emptyset, later fields override earlier fields; giving a struct in
%             place of a name-value pair is equivalent to writing out all the struct fieldnames and
%             respective values. fields that can be derived from others are derived.
%
%   `         optional special semantics:
%             * 'chanlocs' can be specified as cell-string array, and is generally completed using a 
%                default lookup
%             * 'data' can be specified as a cell array of data arrays, then concatenated across 
%                epochs, and with .epoch.target derived as the index of the cell which contained the 
%                epoch in question.
%             * 'tracking.online_expression' can be specified to override the online processing 
%                description
%
% Out:
%   Dataset : newly created EEGLAB set
% 
% Example:
%   % create a new continuous data set (with channels A, B, and C, and 1000 Hz sampling rate)
%   myset = set_new('data',randn(3,100000), 'srate',1000,'chanlocs',struct('labels',{'A','B','C'}));
%
%   % as before, but now also put in some events at certain latencies (note: latencies are in samples)
%   events = struct('type',{'X','Y','X','X','Y'},'latency',{1000,2300,5000,15000,17000});
%   myset = set_new('data',randn(3,100000), 'srate',1000, 'chanlocs',struct('labels',{'A','B','C'}), 'event',events);
%
% See also:
%   eeg_emptyset
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-28

% set_new_version<1.0> -- for the cache

if ~exp_beginfun('filter') return; end

declare_properties('independent_channels',false,'independent_trials',false);

% construct the data set from arguments and defaults
signal = hlp_varargin2struct(varargin, eeg_emptyset, 'setname','new set');

% rewrite cell array data
if iscell(signal.data)
    data = signal.data;
    signal.data = [];
    for c=1:length(data)
        if isnumeric(data{c}) && ~isempty(data{c})
            signal.data = cat(3, signal.data, data{c});
            for i=size(signal.data,3)-size(data{c},3)+1:size(signal.data,3)
                signal.epoch(i).target = c; end
        end
    end    
end

% bring chanlocs into an appropriate format
try signal.chanlocs = hlp_microcache('setnew',@set_infer_chanlocs,signal.chanlocs); catch,end

% if necessary, create chanlocs from scratch, according to the data size
if ~isfield(signal,'chanlocs') || isempty(signal.chanlocs)
    if ~isempty(signal.data)
        signal.chanlocs = struct('labels',cellfun(@num2str,num2cell(1:size(signal.data,1),1),'UniformOutput',false),'type',repmat({'unknown'},1,size(signal.data,1))); 
    else
        signal.chanlocs = struct('labels',{},'type',{});
    end
end

% derive .xmax, .nbchan, .pnts, .trials
[signal.nbchan,signal.pnts, signal.trials, extra_dims] = size(signal.data); %#ok<NASGU>
signal.xmax = signal.xmin + (signal.pnts-1)/signal.srate;

% if epoched and there are events, derive the .epoch field
if signal.trials > 1 && ~isempty(signal.event) && isempty(signal.epoch)
    signal = eeg_checkset(signal,'eventconsistency'); end

% add .epoch.latency if possible
if ~isfield(signal.epoch,'latency')
    for i=1:length(signal.epoch)
        try
            tle = [signal.epoch(i).eventlatency{:}]==0;
            if any(tle)
                signal.epoch(i).latency = b.event(b.epoch(i).event(tle)).latency; end
        catch
        end
    end
end

% create .urevent field if applicable
if isempty(signal.urevent) && ~isempty(signal.event)
    signal.urvent = signal.event;
    [signal.event.urevent] = arraydeal(1:length(signal.event));
end

% do minimal consistency checks
if ~isempty(signal.chanlocs) && ~isempty(signal.data) && (length(signal.chanlocs) ~= signal.nbchan)
    if length(signal.chanlocs) == signal.nbchan+1 
        if isfield(signal,'ref') && isscalar(signal.ref) && signal.ref <= signal.nbchan
            signal.chanlocs = signal.chanlocs(setdiff(1:signal.nbchan,signal.ref));
        elseif any(strcmpi({signal.chanlocs.labels},'ref'))
            signal.chanlocs = signal.chanlocs(~strcmpi({signal.chanlocs.labels},'ref'));
        else
            error('The number of supplied channel locations does not match the number of channels in the data. Apparently includes 1 special channel...'); 
        end
    else
        error('The number of supplied channel locations does not match the number of channels in the data.'); 
    end
end
if isfield(signal,'epoch') && ~isempty(signal.epoch) && length(signal.epoch) ~= size(signal.data,3)
    error('The number of data epochs does not match the number of entries in the epoch field'); end

% ensure that .tracking.timeseries_fields is present
if ~isfield(signal,'tracking') || ~isfield(signal.tracking,'timeseries_fields')
    signal.tracking.timeseries_fields = {}; end

if isfield(signal,'tracking') && isfield(signal.tracking,'online_expression')
    % if an online expression was explicitly assigned in set_new, use that
    exp_endfun('set_online',signal.tracking.online_expression);
else
    % otherwise, we treat this as raw data
    exp_endfun('set_online',struct('head',@rawdata,'parts',{{{signal.chanlocs.labels},unique({signal.chanlocs.type})}}));
end
