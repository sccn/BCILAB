function signal = set_selepos(varargin)
% Select a set of epochs from a data set.
% Signal = set_selepos(Signal, EpochIndices)
%
% In:
%   Signal       : epoched EEGLAB data set from which epochs shall be selected
%
%   EpochIndices : indices of the epochs that should be retained
% 
% Out:
%   Signal       : Newly created data set consisting of the retained epochs in the order
%                  specified in EpochIndices.
%
% Examples:
%   % for an epoched data set, select epochs 50-100
%   eeg = set_selepos(eeg,50:100)
%
% See also:
%   set_selinterval, set_partition, set_makepos
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-04-01
dp;

% set_selepos_version<2.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'epoch_range','EpochIndices','epos'},[],uint32([1 1000000]),'Indices of retained epochs.','shape','row'));

% input validation
utl_check_fields(signal,{'epoch','event','data','pnts'},'signal','signal'); %#ok<NODEF>
if isempty(signal.epoch) && size(signal.data,3) > 1
    error('Your signal has epoched data but has an empty .epoch field. This is not permitted. Use set_makepos instead of pop_epoch for epoching in BCILAB.'); end

if ~isempty(signal.epoch)
    if ~isfield(signal.epoch,'event')
        error('The given signal is lacking the required .epoch.event field.'); end
    if isa(epoch_range,'logical') %#ok<NODEF>
        epoch_range = find(epoch_range); end
    if ~isreal(epoch_range)
        error('The given epoch range contains data that is not real-valued.'); end
    if length(signal.epoch) < max(epoch_range)
        error('The epoch range to select exceeds the length of the signal''s .epoch field (%i): %s.',length(signal.epoch),hlp_tostring(epoch_range)); end

    % select epochs within time series fields
    for field = utl_timeseries_fields(signal) 
        if ~isempty(signal.(field{1}))
            try
                signal.(field{1}) = signal.(field{1})(:,:,epoch_range,:,:,:,:,:); 
            catch e
                error('The given epoch indices could not be extracted from the time-series field .%s with error: %s. Field size was: %i, epoch indices were: %s',field{1},e.message,size(signal.(field{1}),3),hlp_tostring(epoch_range));
            end
        end
    end

    % reindex .epoch field
    signal.epoch = signal.epoch(epoch_range);
    % reindex .event field
    signal.event = signal.event([signal.epoch.event]);
    if ~isempty(epoch_range)
        % get new .event.epoch indices
        epoch_numevents = cellfun('length',{signal.epoch.event});
        for j=length(epoch_numevents):-1:1
            tmp{j} = ones(1,epoch_numevents(j))*j; end;
        epoch_indices = [tmp{:}];
        if ~isempty(signal.event)
            % validate signal.event
            if ~isfield(signal.event,'latency')
                error('The given signal is lacking the required .event.latency field.'); end
            if ~isfield(signal.event,'epoch')
                error('The given signal is lacking the required .event.epoch field.'); end
            latency_numels = cellfun('prodofsize',{signal.event.latency});
            if any(latency_numels == 0)
                error('The given signal has one or more events with empty .latency field. This is not permitted.');
            elseif any(latency_numels ~= 1)
                error('The given signal has one or more events with a .latency value that is not a scalar. This is not permitted.');
            end
            % rewrite pseudo-continuous event latencies based on epoch index change
            [signal.event.latency] = arraydeal([signal.event.latency] + signal.pnts*(epoch_indices - [signal.event.epoch]));
            % relink events and epochs
            [signal.event.epoch] = arraydeal(epoch_indices);
        end
        [signal.epoch.event] = chopdeal(1:length(signal.event),epoch_numevents);
    end
    % update misc meta-data
    signal.trials = size(signal.data,3);
    if ~signal.trials
        disp_once('WARNING: After set_selepos your signal is empty.'); end
elseif ~isempty(epoch_range)
    error('The given epochs could not be selected since the signal is empty.');
end

exp_endfun;
