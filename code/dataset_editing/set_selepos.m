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

% set_selepos_version<2.0> -- for the cache

if ~exp_beginfun('editing') return; end

declare_properties('independent_channels',true,'independent_trials',false);

arg_define(varargin, ...
    arg_norep({'signal','Signal'}), ...
    arg({'epoch_range','EpochIndices','epos'},[],[],'Indices of retained epochs.','shape','row'));

if isa(epoch_range,'logical') %#ok<NODEF>
    epoch_range = find(epoch_range); end

% select epochs within time series fields
for field = utl_timeseries_fields(signal) %#ok<NODEF>
    if ~isempty(signal.(field{1}))
        signal.(field{1}) = signal.(field{1})(:,:,epoch_range,:,:,:,:,:); end
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
        % rewrite pseudo-continuous event latencies based on epoch index change
        [signal.event.latency] = arraydeal([signal.event.latency] + signal.pnts*(epoch_indices - [signal.event.epoch]));
        % relink events and epochs
        [signal.event.epoch] = arraydeal(epoch_indices);
    end
    [signal.epoch.event] = chopdeal(1:length(signal.event),epoch_numevents);
end
% update misc meta-data
signal.trials = size(signal.data,3);

exp_endfun;
